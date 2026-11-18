"""Independent missense RNA evidence; no clinical rescue or scoring rule."""
from dataclasses import dataclass, replace

from .routing import (EvidenceDecision, EvidenceRequest, canonical_bytes,
                      canonical_sha, require_sha)
from .spliceai_source import IndexedSpliceSource, digest
from .execution import require_current_policy

KINDS = ('AG','AL','DG','DL')
FIELDS = tuple(f'{prefix}_{k}' for prefix in ('DS','DP') for k in KINDS) + tuple(f'event_{k}_1based' for k in KINDS)
UNITS = {k: 'delta_score' if k.startswith('DS_') else 'genomic_base_offset' if k.startswith('DP_') else 'genomic_coordinate_1based' for k in FIELDS}


@dataclass(frozen=True)
class SpliceRecipe:
    assembly: str
    reference_sha256: str
    annotation_sha256: str
    source_manifest_sha256: str
    release: str
    model: str
    annotation: str
    masked: bool
    distance: int
    reference_hash_key: str = 'genome_manifest_sha256'
    engine: str = 'precomputed'
    runner_id: str | None = None
    runner_sha256: str | None = None
    schema: str = 'spliceai-recipe/v1'

    def __post_init__(self):
        for x in (self.reference_sha256,self.annotation_sha256,self.source_manifest_sha256): require_sha(x)
        if (self.schema != 'spliceai-recipe/v1' or self.assembly != 'GRCh38'
            or not all(isinstance(v,str) and v for v in (self.release,self.model,self.annotation))
            or type(self.masked) is not bool or type(self.distance) is not int or self.distance <= 0
            or self.engine not in ('precomputed','local')
            or self.reference_hash_key not in ('genome','reference_sha256','chromosome_fasta_sha256','genome_manifest_sha256')):
            raise ValueError('Explicit versioned SpliceAI recipe required')
        if self.engine == 'local':
            require_sha(self.runner_sha256)
            if not self.runner_id: raise ValueError('Local recipe requires a qualified runner identity')
        elif self.runner_id is not None or self.runner_sha256 is not None:
            raise ValueError('Precomputed and local recipes cannot be pooled')

    @classmethod
    def from_manifest(cls, manifest):
        fields=('assembly','reference_sha256','annotation_sha256','reference_hash_key','release','model','annotation','masked','distance')
        return cls(**{k:manifest.data[k] for k in fields},source_manifest_sha256=manifest.sha256)

    @property
    def sha256(self): return canonical_sha(self)


class LocalSpliceRunnerAdapter:
    """Bind an existing QualifiedRunner spec; never install or authorize it.

    The pin is canonical_sha({'schema':'spliceai-local-runner/v1',
    'spec': qualified_runner.spec}); its file pins include weights/annotation.
    Actual source/recipe qualification remains a prerequisite for configuration.
    """
    def __init__(self, qualified_runner, *, expected_sha256):
        require_sha(expected_sha256)
        actual=canonical_sha({'schema':'spliceai-local-runner/v1','spec':qualified_runner.spec})
        if actual != expected_sha256: raise ValueError('Qualified local runner spec pin differs')
        self.runner=qualified_runner
        self.runner_id=qualified_runner.spec.runner_id
        self.runner_sha256=actual

    def run(self,request,*,budget):
        if budget is None: raise ValueError('Explicit invocation execution budget required')
        return self.runner.run(request,budget=budget)


class SpliceEvidenceProvider:
    provider_id = 'rna.spliceai'

    def __init__(self, recipe, source=None, *, runner=None, budget=None):
        self.recipe=recipe; self.source=source; self.runner=runner; self.budget=budget

    def plan(self, context):
        """Every genomic SNV is considered; proximity does not exclude exons."""
        v=context.variant
        identity={'assembly':v.assembly if v else self.recipe.assembly}
        reason=None
        if v is None or not context.genomic_validated:
            reason='unqualified_genomic_identity: '+context.status
        else:
            identity.update(chromosome=v.chromosome,position_1based=v.position_1based,ref=v.ref,alt=v.alt,
                            reference_sha256=self.recipe.reference_sha256,annotation_sha256=self.recipe.annotation_sha256)
            if v.assembly != self.recipe.assembly: reason='unsupported_assembly'
            elif len(v.ref) != 1 or len(v.alt) != 1 or v.ref not in 'ACGT' or v.alt not in 'ACGT' or v.ref == v.alt:
                reason='unsupported_non_snv'
            elif v.position_1based != v.position_0based+1: reason='inconsistent_normalized_coordinate'
        transcripts={t.transcript_id:t for t in context.transcripts}
        selected=tuple(sorted({m.transcript_id for m in context.mappings}))
        ids=selected or tuple(sorted(transcripts)) or (None,)
        result=[]
        for tx in ids:
            item=dict(identity); params={'route_reason':reason or 'eligible_genomic_snv',
                'applicable':reason is None,'source_manifest_sha256':self.recipe.source_manifest_sha256,
                'coordinate_normalization':'normalized_position_0based_plus_one',
                'proximity_filter':'none'}
            if tx is not None:
                item['transcript_id']=tx
                t=transcripts.get(tx)
                m=next((m for m in context.mappings if m.transcript_id==tx),None)
                strand=t.strand if t else m.strand if m else None
                if strand: item['strand']=strand
                if t is not None:
                    params.update(reference_hashes=dict(t.reference_hashes),boundary_status=t.boundary_status)
                    # These are context only; no prediction is excluded by distance.
                    params['exon_intervals']=t.exon_intervals
                else: params['reference_hashes']={}
            result.append(EvidenceRequest(self.provider_id,self.recipe.sha256,item,params))
        return tuple(result)

    def _decision(self, request, status, reason, **kw):
        return EvidenceDecision(request,status != 'not_applicable',status,reason,**kw)

    def _source_contract(self):
        if self.source is None: return
        m=self.source.manifest
        if m.sha256 != self.recipe.source_manifest_sha256:
            raise ValueError('Source manifest pin differs from recipe')
        for field in ('assembly','reference_sha256','annotation_sha256','reference_hash_key','release','model','annotation','masked','distance'):
            if m.data[field] != getattr(self.recipe,field): raise ValueError('Source recipe conflict: '+field)

    def _compatibility(self,request):
        tx=request.identity.get('transcript_id')
        if not tx: return None,'selected_transcript_unqualified'
        refs=request.parameters.get('reference_hashes',{})
        genome=refs.get(self.recipe.reference_hash_key)
        if genome != self.recipe.reference_sha256:
            return None,'reference_hash_unqualified_or_conflicting'
        if 'annotation' in refs and refs['annotation'] != self.recipe.annotation_sha256:
            return None,'annotation_hash_conflict'
        if self.source is None: return None,'source_absent'
        rows=[r for r in self.source.manifest.data['mane_records'] if r['transcript_id']==tx]
        if len(rows) != 1: return None,'versioned_mane_transcript_unqualified'
        mane=rows[0]
        if request.identity.get('strand') != mane['strand']: return None,'mane_strand_conflict'
        if request.identity.get('gene_id',mane['gene_id']) != mane['gene_id'] or request.identity.get('gene_symbol',mane['gene_symbol']) != mane['gene_symbol']:
            return None,'mane_gene_conflict'
        return mane,'exact_versioned_mane'

    def obtain(self, request, store, policy):
        if request.provider_id != self.provider_id or request.recipe_sha256 != self.recipe.sha256:
            return self._decision(request,'failed','provider_or_recipe_pin_mismatch')
        if not request.parameters.get('applicable',False):
            return self._decision(request,'not_applicable',request.parameters.get('route_reason','unqualified_route'))
        # Validate source/request identity before trusting an exact cached record.
        try:
            i=request.identity
            for field in ('assembly','reference_sha256','annotation_sha256'):
                if i.get(field) != getattr(self.recipe,field):
                    raise ValueError('Request identity differs from pinned recipe: '+field)
            if request.parameters.get('source_manifest_sha256') != self.recipe.source_manifest_sha256:
                raise ValueError('Request source manifest differs from recipe')
            if (type(i.get('position_1based')) is not int or i['position_1based'] < 1
                or i.get('ref') not in ('A','C','G','T') or i.get('alt') not in ('A','C','G','T')
                or i['ref'] == i['alt']): raise ValueError('Request lacks an exact genomic SNV')
            self._source_contract()
        except ValueError as exc: return self._decision(request,'failed',str(exc))
        if store is not None:
            cached=store.lookup(request)
            if cached is not None: return cached
        def produce():
            if self.recipe.engine == 'local': return self._run_local(request,policy)
            return self._obtain_source(request,policy)
        if store is not None and hasattr(store,'get_or_compute'):
            return store.get_or_compute(request,produce)
        artifact,result=produce()
        if store is not None and result.status in ('available','partial','ambiguous','unavailable') and result.metadata.get('source_queried'):
            return store.publish(request,artifact,result)
        return result

    def _obtain_source(self,request,policy):
        empty=b''
        if self.source is None: return empty,self._decision(request,'unavailable','qualified_score_source_not_configured')
        mane,compat=self._compatibility(request)
        if mane is None:
            return empty,self._decision(request,'ambiguous',compat,metadata={'source_scope':'gene','transcript_compatibility':compat})
        remote=isinstance(self.source,IndexedSpliceSource) or getattr(self.source,'network_required',False)
        if remote and policy.mode != 'bounded-retrieval':
            return empty,self._decision(request,'unavailable','cache_only_no_qualified_cached_score',cache_disposition='miss')
        try:
            if remote:require_current_policy(policy,fetcher=getattr(self.source,'fetcher',None))
            i=request.identity
            source_result=self.source.query(i['chromosome'],i['position_1based'],i['ref'],i['alt'],variant_key=canonical_sha({k:v for k,v in i.items() if k in ('assembly','chromosome','position_1based','ref','alt')}))
            candidates=[r for r in source_result.records if r.gene_symbol==mane['gene_symbol']]
            metadata={'source_scope':'gene','transcript_compatibility':compat,'gene_id':mane['gene_id'],
                'gene_symbol':mane['gene_symbol'],'selected_transcript_id':i['transcript_id'],
                'annotation':self.recipe.annotation,'source_queried':True,'event_coordinate_rule':'VCF_POS_plus_DP_no_strand_inversion',
                'coverage':'exact_allele_gene','range_receipts':source_result.receipts}
            provenance={'source_manifest_sha256':self.recipe.source_manifest_sha256,
                'recipe_sha256':self.recipe.sha256,'source_payload_sha256':digest(source_result.artifact_bytes),
                'reference_sha256':self.recipe.reference_sha256,'annotation_sha256':self.recipe.annotation_sha256}
            if len(candidates) != 1:
                status='ambiguous' if candidates else 'unavailable'
                reason='ambiguous_gene_score' if candidates else 'missing_exact_allele_gene_score'
                return source_result.artifact_bytes,self._decision(request,status,reason,metadata=metadata,provenance=provenance)
            record=candidates[0]
            if any(v is not None and abs(v) > self.recipe.distance for v in record.offsets):
                raise ValueError('Source offset exceeds pinned recipe distance')
            values={f'DS_{k}':s for k,s in zip(KINDS,record.scores)}
            values.update({f'DP_{k}':s for k,s in zip(KINDS,record.offsets)})
            values.update({f'event_{k}_1based':None if o is None else i['position_1based']+o for k,o in zip(KINDS,record.offsets)})
            if any(v is not None and v < 1 for k,v in values.items() if k.startswith('event_')):
                raise ValueError('Predicted event coordinate outside reference support')
            known=sum(v is not None for v in values.values())
            status='available' if known == len(values) else 'partial' if known else 'unavailable'
            return source_result.artifact_bytes,self._decision(request,status,'qualified_gene_scope_prediction' if known else 'all_source_scores_missing',
                values=values,units=UNITS,provenance=provenance,metadata=metadata,
                uncertainty='predicted_splicing_not_measured_RNA_or_transcript_specific_effect',
                execution='bounded_retrieval' if remote else 'local_artifact_read')
        except (ValueError,OSError,RuntimeError) as exc:
            return empty,self._decision(request,'failed','source_intake_failed: '+str(exc),execution='source_read_failed')

    def _run_local(self,request,policy):
        if self.runner is None:
            return b'',self._decision(request,'unavailable','qualified_local_spliceai_runner_or_weights_or_annotation_absent')
        refs=request.parameters.get('reference_hashes',{})
        if refs.get(self.recipe.reference_hash_key)!=self.recipe.reference_sha256:
            return b'',self._decision(request,'ambiguous','reference_hash_unqualified_or_conflicting')
        if 'annotation' in refs and refs['annotation']!=self.recipe.annotation_sha256:
            return b'',self._decision(request,'ambiguous','annotation_hash_conflict')
        if policy.mode != 'bounded-retrieval' or self.recipe.runner_id not in policy.allowed_runner_ids or not policy.authorization_refs:
            return b'',self._decision(request,'pending_authorization','current_local_runner_authorization_required')
        if (getattr(self.runner,'runner_id',None),getattr(self.runner,'runner_sha256',None)) != (self.recipe.runner_id,self.recipe.runner_sha256):
            return b'',self._decision(request,'failed','local_runner_pin_mismatch')
        try:
            require_current_policy(policy,budget=self.budget)
            generic=getattr(self.runner,'runner',None)
            if generic is not None and getattr(generic,'policy',policy)!=policy:
                raise ValueError('Configured local runner policy differs from current policy')
            result=self.runner.run(request,budget=self.budget)
            if not isinstance(result,EvidenceDecision) or result.request != request:
                raise ValueError('Runner returned incompatible request')
            if result.status in ('available','partial') and (set(result.values) != set(FIELDS) or result.scoring_used):
                raise ValueError('Runner must return four scores, offsets and genomic events without scoring')
            if result.status in ('available','partial'):
                for k in KINDS:
                    score,offset,event=(result.values[f'DS_{k}'],result.values[f'DP_{k}'],result.values[f'event_{k}_1based'])
                    if score is not None and not 0 <= score <= 1: raise ValueError('Local delta score outside 0..1')
                    if offset is not None and (type(offset) is not int or abs(offset) > self.recipe.distance):
                        raise ValueError('Local offset outside pinned recipe')
                    expected=None if offset is None else request.identity['position_1based']+offset
                    if event != expected or event is not None and (type(event) is not int or event < 1):
                        raise ValueError('Local event does not equal VCF POS plus DP')
                if dict(result.units) != UNITS: raise ValueError('Local runner units differ from recipe')
                result=replace(result,qualification='diagnostic',scoring_used=(),
                    metadata=dict(result.metadata,source_scope='gene',engine='local',source_queried=True),
                    provenance=dict(result.provenance,recipe_sha256=self.recipe.sha256,runner_sha256=self.recipe.runner_sha256))
            return canonical_bytes(result),result
        except (ValueError,OSError,RuntimeError) as exc:
            return b'',self._decision(request,'failed','local_runner_failed: '+str(exc))
