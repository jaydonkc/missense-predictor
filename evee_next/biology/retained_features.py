"""Complete retained-coordinate outputs and fixed training-transform application.

This module neither extracts backbone features nor qualifies a recipe by itself.
An explicitly pinned compatibility receipt is required for supplied outputs.
Legacy query vectors are not a representation of the retained model.
"""
from io import BytesIO
from pathlib import Path
import hashlib,zipfile
import numpy as np
from evee_next.report import MissenseMapping,NormalizedVariant
from evee_next.interpretation.neural import immutable
from .admission import pinned_json
from .neural_pipeline import NeuralFeatureRecord
from .retained_provider import BLOCKS,SOURCES
from .routing import require_sha,canonical_sha


def apply_frozen_transform(raw,mean,scale,projection):
    """Apply supplied fitting-only statistics; never infer statistics from queries."""
    x,m,s,p=(np.asarray(v,dtype=np.float64) for v in (raw,mean,scale,projection))
    if x.ndim!=1 or m.shape!=x.shape or s.shape!=x.shape or p.ndim!=2 or p.shape[0]!=len(x):
        raise ValueError('Fixed transformation dimensions differ')
    if not all(np.isfinite(v).all() for v in (x,m,s,p)) or (s<=0).any():
        raise ValueError('Finite fitting-only statistics and positive scales required')
    result=((x-m)/s)@p
    if not np.isfinite(result).all():raise ValueError('Nonfinite transformed features')
    return immutable(result)


class RetainedFeatureProvider:
    def __init__(self,manifest_path,manifest_sha256,*,checkpoint):
        root=Path(manifest_path).parent;d=pinned_json(manifest_path,manifest_sha256)
        if set(d)!={'schema','stage','checkpoint','feature_map_id','qualification','artifact','recipe_hashes','rows'} or d['schema']!='retained-feature-outputs/v1' or d['stage']!='original_fitting_only_transformed_coordinates':
            raise ValueError('Complete retained-coordinate output contract required')
        if (d['checkpoint'],d['feature_map_id'])!=(checkpoint.parameter_sha256,checkpoint.feature_map_id):
            raise ValueError('Output checkpoint/coordinate map differs')
        if tuple((b.name,b.start,b.end) for b in checkpoint.blocks)!=tuple((b.name,b.start,b.end) for b in BLOCKS):
            raise ValueError('Complete original eight-block retained model required')
        if set(d['recipe_hashes'])!={b.name for b in BLOCKS}:raise ValueError('Every original block recipe must be pinned')
        for v in d['recipe_hashes'].values():require_sha(v)
        q=d['qualification'];receipt=pinned_json(root/q['path'],q['sha256'])
        if (receipt.get('schema'),receipt.get('status'),receipt.get('kind'),receipt.get('checkpoint'),receipt.get('feature_map_id')) != (
            'retained-feature-qualification/v1','passed','complete_feature_compatibility',d['checkpoint'],d['feature_map_id']):
            raise ValueError('Passed checkpoint-compatible extraction receipt required')
        artifact=d['artifact'];require_sha(artifact['sha256'])
        if (receipt.get('artifact_sha256'),receipt.get('recipe_hashes_sha256'),receipt.get('rows_sha256')) != (
            artifact['sha256'],canonical_sha(d['recipe_hashes']),canonical_sha(d['rows'])):
            raise ValueError('Feature qualification does not bind exact recipes, output bytes and identity roster')
        data=(root/artifact['path']).read_bytes()
        if hashlib.sha256(data).hexdigest()!=artifact['sha256']:raise ValueError('Transformed feature bytes changed')
        with zipfile.ZipFile(BytesIO(data)) as z:
            if z.namelist()!=['X.npy']:raise ValueError('Only complete transformed X belongs in inference output')
        with np.load(BytesIO(data),allow_pickle=False) as z:self.X=immutable(z['X'])
        if self.X.shape!=(len(d['rows']),5333):raise ValueError('Expected complete retained 5,333-column inputs')
        self.rows={}
        self.references={}
        for i,row in enumerate(d['rows']):
            if set(row)!={'variant','mapping','reference_hashes'}:raise ValueError('Outcome-free exact output identity and reference hashes required')
            refs=row['reference_hashes']
            if not isinstance(refs,dict) or 'protein_sha256' not in refs or not {'genome_manifest_sha256','chromosome_fasta_sha256'}.intersection(refs):
                raise ValueError('Full protein and genomic source hashes required')
            for value in refs.values():require_sha(value)
            v,m=NormalizedVariant(**row['variant']),MissenseMapping(**row['mapping'])
            from evee_next.report import VariantRequest,normalize
            if normalize(VariantRequest(v.chromosome,v.position_1based,v.ref,v.alt,assembly=v.assembly,organism=v.organism))!=v:
                raise ValueError('Feature output normalized identity differs')
            key=(v,m)
            if key in self.rows:raise ValueError('Duplicate transformed output identity')
            self.rows[key]=i
            self.references[key]=dict(refs)
        self.feature_map_id=d['feature_map_id']
        self.provenance={'feature_output_manifest_sha256':manifest_sha256,
            'feature_qualification_sha256':q['sha256'],'feature_scope':'explicit_complete_output_compatibility_not_general_extraction'}

    def validate_context(self,context):
        key=(context.variant,context.mapping)
        if key not in self.references:return False
        expected=self.references[key]
        actual=dict(context.provenance)
        for t in context.transcripts:
            if (t.transcript_id,t.protein_id)==(context.mapping.transcript_id,context.mapping.protein_id):
                actual.update(t.reference_hashes)
        from .baseline import context_protein_sha256
        actual['protein_sha256']=context_protein_sha256(context)
        if any(actual.get(k)!=v for k,v in expected.items()):
            raise ValueError('Current context reference/source hashes differ from qualified feature outputs')
        return True

    def features(self,variant,mapping,required_blocks):
        if tuple(required_blocks)!=tuple(b.name for b in BLOCKS):raise ValueError('Original retained block order required')
        i=self.rows.get((variant,mapping))
        values={} if i is None else {b.name:self.X[i,b.start:b.end] for b in BLOCKS}
        return NeuralFeatureRecord(variant,mapping,self.feature_map_id,values,SOURCES,self.provenance)
