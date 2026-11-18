"""Checkpoint-bound retained inference, preserving original cached fold ownership."""
from pathlib import Path
from evee_next.report import normalize, matching_mappings
from evee_next.reference import TranscriptContext
from .retained_provider import CachedRetainedProvider, _manifest, SCOPE
from .routing import BaselinePrediction,BaselineResult,BaselineState,VariantContext


def context_protein_sha256(context):
    import hashlib
    m=context.mapping
    rows=[t for t in context.transcripts if t.transcript_id==m.transcript_id and t.protein_id==m.protein_id]
    if len(rows)!=1:raise ValueError('Unique full protein reference context required')
    row=rows[0];claimed=row.reference_hashes.get('protein_sha256')
    if row.protein_sequence is not None:
        actual=hashlib.sha256(row.protein_sequence.encode()).hexdigest()
        if claimed is not None and claimed!=actual:raise ValueError('Context protein sequence/hash conflict')
        return actual
    if claimed is None:raise ValueError('Full protein reference hash unavailable')
    return claimed


class CachedBaselineProvider:
    def __init__(self,bundle,manifest_sha256):
        self.bundle=Path(bundle);self.pin=manifest_sha256
        self.manifest=_manifest(self.bundle,self.pin)
        self.providers={}
        self._folds={}
        for row in self.manifest['rows']:
            key=row['variant']['variant_id'];fold=row['fold']
            if self._folds.setdefault(key,fold)!=fold:
                raise ValueError('An archived allele cannot borrow another original head')

    def _provider(self,fold):
        if fold not in self.providers:
            self.providers[fold]=CachedRetainedProvider(self.bundle,self.pin,fold)
        return self.providers[fold]

    def fold_for(self,variant): return self._folds.get(variant.variant_id)

    def context_for(self,request):
        try: variant=normalize(request)
        except NotImplementedError as exc:
            return VariantContext('unsupported',str(exc),None)
        except ValueError as exc:
            return VariantContext('invalid_request',str(exc),None)
        fold=self.fold_for(variant)
        if fold is None:
            return VariantContext('outside_admitted_roster','No archived mapping; provide a qualified fresh resolver for novel diagnostic queries',variant)
        provider=self._provider(fold)
        options=provider.mappings(variant);chosen=matching_mappings(request,options)
        if not chosen:
            return VariantContext('mapping_conflict','Requested version does not match archived identity',variant,
                                  options,provenance=provider.provenance)
        contexts=[]
        for m in chosen:
            ref=provider.references[m.protein_id]
            contexts.append(TranscriptContext(m.transcript_id,m.protein_id,m.strand,(),(),
                {'protein_sha256':ref['sha256'],'retained_manifest_sha256':self.pin},'unavailable',ref['sequence']))
        return VariantContext('resolved' if len(chosen)==1 else 'ambiguous_mapping',
            'Archived qualified mapping; no fresh genomic REF/codon remap',variant,chosen,
            chosen[0] if len(chosen)==1 else None,tuple(contexts),provider.provenance,True)

    def predict(self,context):
        if context.variant is None or context.mapping is None:
            return BaselineResult(BaselinePrediction('unscored',reason='Unique qualified mapping required'))
        fold=self.fold_for(context.variant)
        if fold is None:
            return BaselineResult(BaselinePrediction('outside_admitted_roster',reason='No complete retained feature coordinates for this allele'))
        provider=self._provider(fold);head=provider.head
        if context_protein_sha256(context)!=provider.references[context.mapping.protein_id]['sha256']:
            raise ValueError('Current full protein reference differs from retained cached coordinates')
        record=provider.features(context.variant,context.mapping,tuple(b.name for b in head.blocks))
        if set(record.blocks)!={b.name for b in head.blocks}:
            return BaselineResult(BaselinePrediction('partial',reason='Complete original feature blocks required'))
        import numpy as np
        x=np.concatenate([record.blocks[b.name] for b in head.blocks])
        score=head.account(x,record.feature_map_id)
        prediction=BaselinePrediction('scored',score.value,head.cutoff,head.parameter_sha256,
            head.feature_map_id,SCOPE,reason='retained_development; original held-out head; uncalibrated score')
        state=BaselineState(context.variant.variant_id,head.parameter_sha256,head.feature_map_id,
            tuple(float(v) for v in (x@head.W0+head.b0)),tuple(float(v) for v in head.W1),head.b1)
        return BaselineResult(prediction,state)


class QualifiedCheckpointBaseline:
    """Explicit development checkpoint for complete qualified query outputs.

    Caller must supply a fresh identity resolver separately. This does not pick
    an assessment fold, refit a head, ensemble maps or extract raw sequences.
    """
    def __init__(self,head,feature_provider):
        self.head,self.feature_provider=head,feature_provider

    def predict(self,context):
        import numpy as np
        if not context.genomic_validated or context.mapping is None:
            return BaselineResult(BaselinePrediction('unscored',reason='Unique qualified reference identity required'))
        head=self.head
        if not hasattr(self.feature_provider,'validate_context'):
            raise ValueError('Qualified outputs need a current reference/source binding')
        if not self.feature_provider.validate_context(context):
            return BaselineResult(BaselinePrediction('outside_qualified_feature_scope',reason='No identity/reference-bound output for this allele'))
        record=self.feature_provider.features(context.variant,context.mapping,tuple(b.name for b in head.blocks))
        if (record.variant,record.mapping,record.feature_map_id)!=(context.variant,context.mapping,head.feature_map_id):
            raise ValueError('Qualified output identity/coordinate map differs')
        if set(record.blocks)!={b.name for b in head.blocks}:
            return BaselineResult(BaselinePrediction('outside_qualified_feature_scope',reason='No complete checkpoint-compatible output for this allele'))
        x=np.concatenate([record.blocks[b.name] for b in head.blocks]);score=head.account(x,record.feature_map_id)
        prediction=BaselinePrediction('scored',score.value,head.cutoff,head.parameter_sha256,head.feature_map_id,
            'explicit_checkpoint_development',reason='qualified supplied feature outputs; not a serving fit or independently confirmed predictor')
        return BaselineResult(prediction,BaselineState(context.variant.variant_id,head.parameter_sha256,head.feature_map_id,
            tuple(x@head.W0+head.b0),tuple(head.W1),head.b1))


class CachedThenQualifiedBaseline:
    """Preserve cached assessment ownership; explicit outputs extend coverage."""
    def __init__(self,cached,qualified):self.cached,self.qualified=cached,qualified
    def context_for(self,request):return self.cached.context_for(request)
    def predict(self,context):
        if self.cached.fold_for(context.variant) is not None:return self.cached.predict(context)
        return self.qualified.predict(context)
