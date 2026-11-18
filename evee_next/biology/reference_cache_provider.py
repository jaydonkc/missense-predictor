"""Reference-only projection of the existing qualified target-free biology cache.

Reads metadata and context.npy only. Effect heads, protein features, chemistry
arrays and scientific outcomes are never decoded or imported by this adapter.
Precomputed RSA retains its historical source recipe; it is not recomputed or
relabeled as Shrake–Rupley or a mutation effect.
"""
from io import BytesIO
from pathlib import Path
from types import MappingProxyType
import hashlib,json,zipfile
import numpy as np
from .routing import EvidenceDecision,canonical_sha,canonical_bytes,require_sha
from .structure_provider import _request


PROTEINS={'ASPA':'ENSP00000263080.2','G6PD':'ENSP00000377192.3'}


def snapshot(path,pin,limit=256*1024**2):
    require_sha(pin)
    with Path(path).open('rb') as f:data=f.read(limit+1)
    if len(data)>limit or hashlib.sha256(data).hexdigest()!=pin:
        raise ValueError('Qualified context cache source hash/size differs')
    return data


class CachedReferenceContextProvider:
    provider_id='structure.cached'

    def __init__(self,bundle,manifest_sha256):
        self.bundle=Path(bundle);self.pin=manifest_sha256
        manifest=json.loads(snapshot(self.bundle/'manifest.json',self.pin,1024**2))
        if manifest.get('schema')!='biology-target-free-inputs/v1' or set(manifest.get('files',{}))!={'features.npz','heads.json','metadata.json'}:
            raise ValueError('Existing target-free cache contract required')
        meta=json.loads(snapshot(self.bundle/'metadata.json',manifest['files']['metadata.json'],16*1024**2))
        if meta.get('targets_present') is not False or not isinstance(meta.get('biological_development_lead'),bool):
            raise ValueError('Explicit target-free source qualification required')
        rows,refs=meta.get('rows'),meta.get('references')
        if not isinstance(rows,list) or not rows or not isinstance(refs,dict) or not refs or not set(refs)<=set(PROTEINS):
            raise ValueError('Exact supported reference/row roster required')
        for gene,ref in refs.items():
            if not isinstance(ref.get('sequence'),str) or hashlib.sha256(ref['sequence'].encode()).hexdigest()!=ref.get('sha256'):
                raise ValueError('Full cached reference sequence/hash differs')
            if ref.get('protein_id',PROTEINS[gene])!=PROTEINS[gene]:raise ValueError('Versioned cached protein differs')
        raw=snapshot(self.bundle/'features.npz',manifest['files']['features.npz'])
        with zipfile.ZipFile(BytesIO(raw)) as z:
            if len(z.namelist())!=3 or set(z.namelist())!={'context.npy','chemistry.npy','protein.npy'}:
                raise ValueError('Only the target-free inference archive is supported')
            with z.open('context.npy') as f:context_bytes=f.read(len(rows)*3*8+1024)
        context=np.load(BytesIO(context_bytes),allow_pickle=False)
        if context.dtype.kind!='f' or context.shape!=(len(rows),3) or not np.isfinite(context).all() or ((context<0)|(context>1)).any():
            raise ValueError('Finite normalized pLDDT/RSA/availability context required')
        if not np.isin(context[:,2],[0.,1.]).all() or (context[context[:,2]==0,1]!=0).any():
            raise ValueError('Invalid source availability encoding')
        lookup={};ids=set()
        for row,values in zip(rows,context):
            gene=row.get('gene');position=row.get('position0');rid=row.get('id')
            if gene not in refs or type(position) is not int or not isinstance(rid,str) or not rid or rid in ids:
                raise ValueError('Invalid reference context identity roster')
            ids.add(rid);reference=refs[gene]
            if not 0<=position<len(reference['sequence']) or reference['sequence'][position]!=row.get('ref_aa'):
                raise ValueError('Focal residue differs from complete cached reference')
            if row.get('protein_id',PROTEINS[gene])!=PROTEINS[gene] or row.get('reference_sha256',reference['sha256'])!=reference['sha256']:
                raise ValueError('Row protein/reference binding differs')
            key=(PROTEINS[gene],reference['sha256'],position+1,row['ref_aa'])
            value=tuple(float(x) for x in values)
            if key in lookup and lookup[key]!=value:raise ValueError('Reference context changes across alternate substitutions')
            lookup[key]=value
        self.lookup=MappingProxyType(lookup)
        self.source_sha256=manifest['files']['features.npz']
        self.metadata_sha256=manifest['files']['metadata.json']
        self.recipe_sha256=canonical_sha({'provider':self.provider_id,'manifest_sha256':self.pin,
            'protein_map':PROTEINS,'context_recipe':'original_normalized_plddt_RSA_availability',
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})

    def plan(self,context):return (_request(self.provider_id,self.recipe_sha256,context),)

    def obtain(self,request,store,policy):
        if request.provider_id!=self.provider_id or request.recipe_sha256!=self.recipe_sha256:
            return EvidenceDecision(request,True,'failed','Reference cache provider/recipe differs')
        def produce():
            identity=request.identity
            key=(identity.get('protein_id'),identity.get('protein_sha256'),identity.get('residue_1based'),identity.get('ref_aa'))
            values=self.lookup.get(key)
            if request.parameters.get('identity_status')!='qualified' or values is None:
                return b'',EvidenceDecision(request,True,'unavailable','No exact qualified full-reference/site context in the existing cache')
            q,rsa,available=values
            decision=EvidenceDecision(request,True,'available' if available else 'partial',
                'qualified_precomputed_reference_context',{'pLDDT':q*100,'RSA':rsa if available else None},
                {'pLDDT':'confidence_0_100','RSA':'relative_solvent_accessibility'},
                {'source_manifest_sha256':self.pin,'source_artifact_sha256':self.source_sha256,
                 'metadata_sha256':self.metadata_sha256,'context_recipe':'original_normalized_plddt_RSA_availability'},
                {'reference_context_only':True,'RSA_recipe':'precomputed_source_recipe_not_rederived',
                 'effect_heads_loaded':False,'assay_targets_loaded':False,'mutation_effect_measured':False},
                uncertainty='Reference context; source recipe preserved; not a mutation effect',execution='local_context_projection')
            return canonical_bytes(decision),decision
        if store is not None:
            try:return store.get_or_compute(request,produce)
            except (ValueError,OSError) as exc:
                return EvidenceDecision(request,True,'failed','Reference context cache failed integrity/access: '+str(exc))
        return produce()[1]
