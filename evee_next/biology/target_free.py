"""Portable cached effect diagnostics; no training targets or source workspace.

Requires an explicit manifest SHA trust anchor. The bundle contains features,
references, model parameters and qualification only, never clinical/assay targets.
Coverage is cached ASPA/G6PD substitutions, not general sequence inference.
"""
import argparse
import hashlib
from io import BytesIO
import json
from pathlib import Path
import numpy as np
from .chemistry import Chemistry
from .effects import designs, finite, predict_ridge
from .runtime import PROTEINS, protein_report
from .chemistry import ProteinVariant, TABLE_SHA


def _bound_payload(path,sha256):
    # Decode exactly the bytes that passed the trust anchor, rather than reopening
    # a path after hashing it (the file may have been replaced in between).
    try:payload=path.read_bytes()
    except OSError as error:raise ValueError('Bundle input unavailable: '+str(path)) from error
    if hashlib.sha256(payload).hexdigest()!=sha256:
        raise ValueError('Input hash differs: '+str(path))
    return payload


class TargetFreeEffectProvider:
    def __init__(self,bundle: Path,manifest_sha256: str):
        self.bundle=Path(bundle)
        manifest=json.loads(_bound_payload(self.bundle/'manifest.json',manifest_sha256))
        if manifest['schema']!='biology-target-free-inputs/v1' or set(manifest['files'])!={'features.npz','heads.json','metadata.json'}:
            raise ValueError('Unexpected target-free bundle schema')
        payloads={name:_bound_payload(self.bundle/name,sha) for name,sha in manifest['files'].items()}
        self.meta=json.loads(payloads['metadata.json'])
        self.heads=json.loads(payloads['heads.json'])
        if type(self.meta.get('biological_development_lead')) is not bool or self.meta.get('targets_present',False) is not False:
            raise ValueError('Explicit development qualification and target-free inputs required')
        rows=self.meta['rows'];references=self.meta['references']
        if not isinstance(rows,list) or not rows or not references or not set(references)<=set(PROTEINS):
            raise ValueError('Supported reference and nonempty row roster required')
        self.chemistry=None
        if 'chemistry_tables' in self.meta:
            self.chemistry=Chemistry(self.meta['chemistry_tables'],'pinned_classical_tables#sha256='+TABLE_SHA)
        self.lookup={}
        ids=set()
        for i,row in enumerate(rows):
            gene=row['gene']
            if gene not in PROTEINS or gene not in references:raise ValueError('Unsupported inference source')
            if type(row['position0']) is not int or not isinstance(row['id'],str) or not row['id'] or row['id'] in ids:
                raise ValueError('Unique row ID and exact integer protein position required')
            ids.add(row['id'])
            reference=references[gene]
            if (row.get('protein_id',PROTEINS[gene])!=PROTEINS[gene]
                    or reference.get('protein_id',PROTEINS[gene])!=PROTEINS[gene]
                    or row.get('reference_sha256',reference['sha256'])!=reference['sha256']):
                raise ValueError('Row protein/reference binding differs')
            v=ProteinVariant(gene,PROTEINS[gene],row['position0']+1,row['ref_aa'],row['alt_aa'],reference['sha256'])
            v.validate(reference['sequence'])
            key=(gene,v.position_1based,v.ref_aa,v.alt_aa)
            if key in self.lookup:raise ValueError('Duplicate inference substitution')
            self.lookup[key]=i
        with np.load(BytesIO(payloads['features.npz']),allow_pickle=False) as z:
            if len(z.files)!=3 or set(z.files)!={'chemistry','context','protein'}:
                raise ValueError('Unexpected or duplicate feature members or outcomes in inference bundle')
            c,r,p=(finite(z[name],2) for name in ('chemistry','context','protein'))
        if c.shape!=(len(rows),56) or r.shape!=(len(rows),3) or p.shape!=(len(rows),516):
            raise ValueError('Feature dimensions differ from the released biology recipe')
        if ((r<0)|(r>1)).any() or not np.isin(r[:,2],[0.,1.]).all() or (r[r[:,2]==0,1]!=0).any():
            raise ValueError('Invalid normalized reference context or RSA availability')
        self.designs,_,_=designs(c,r,p,rows)
        self.head_lookup={}
        if not isinstance(self.heads,list):raise ValueError('Saved effect heads must be a list')
        for head in self.heads:
            gene=head['assessment_gene'];arm=head['arm'];key=(gene,arm)
            if gene not in references or arm not in self.designs or key in self.head_lookup:
                raise ValueError('Unexpected or duplicate effect head identity')
            width=self.designs[arm].shape[1]
            fitting_gene='G6PD' if gene=='ASPA' else 'ASPA'
            if (head.get('fitting_gene',fitting_gene)!=fitting_gene
                    or type(head.get('columns',width)) is not int or head.get('columns',width)!=width):
                raise ValueError('Saved head fitting source or feature width differs')
            state=head['state']
            params=[finite(state[name],1) for name in ('mean','scale','coef')]
            if (any(x.shape!=(width,) for x in params) or (params[1]<=0).any()
                    or type(state['intercept']) not in (int,float) or not np.isfinite(state['intercept'])):
                raise ValueError('Invalid saved effect head dimensions, scale or scalar intercept')
            self.head_lookup[key]=head
        for gene in {row['gene'] for row in rows}:
            if any((gene,arm) not in self.head_lookup for arm in ('C','R','K','P','J')):
                raise ValueError('Missing required saved effect head')
        self.context=r;self.manifest_sha256=manifest_sha256

    def __call__(self,v):
        if PROTEINS.get(v.gene)!=v.protein_id:raise ValueError('Exact versioned protein identity required')
        if v.gene not in self.meta['references']:raise ValueError('No admitted reference')
        reference=self.meta['references'][v.gene];v.validate(reference['sequence'])
        if reference['sha256']!=v.reference_sha256:raise ValueError('Reference binding differs')
        i=self.lookup.get((v.gene,v.position_1based,v.ref_aa,v.alt_aa))
        if i is None:return None
        predictions={}
        for arm in ('C','R','K','P','J'):
            head=self.head_lookup[v.gene,arm]
            predictions[arm]=float(predict_ridge(head['state'],self.designs[arm][i:i+1])[0])
        ctx=self.context[i];rsa=bool(ctx[2])
        return dict(endpoint='reporter_abundance',units='fitting_assay_native_fluorescence_bin_scale',
                    fitting_gene='G6PD' if v.gene=='ASPA' else 'ASPA',predictions=predictions,
                    context=dict(reference_pLDDT=float(ctx[0]*100),reference_RSA=float(ctx[1]) if rsa else None,
                                 RSA_available=rsa,source='pinned_reference_context_cache'),
                    qualification='development_lead' if self.meta['biological_development_lead'] else 'failed_advancement_do_not_use_for_clinical_scoring',
                    uncertainty='Not individually calibrated; see held-protein evaluation',
                    source=str((self.bundle/'manifest.json').resolve()),result_sha256=self.manifest_sha256,
                    is_pathogenicity_probability=False)

    def query(self,gene,position,ref,alt):
        if gene not in PROTEINS or gene not in self.meta['references']:raise ValueError('No admitted reference')
        r=self.meta['references'][gene];v=ProteinVariant(gene,PROTEINS[gene],position,ref,alt,r['sha256'])
        chemistry=self.chemistry
        if chemistry is None:
            # Test fixtures can supply an explicit chemistry calculator; real bundles
            # always contain the pinned tables. Never invent chemistry for a query.
            class Missing:
                def calculate(self,ref,alt):return {'status':'unavailable','is_pathogenicity':False}
            chemistry=Missing()
        return protein_report(v,r['sequence'],chemistry,self)


def main():
    p=argparse.ArgumentParser(description='Target-free cached molecular-effect report; no clinical prediction')
    p.add_argument('bundle',type=Path);p.add_argument('manifest_sha256')
    p.add_argument('gene',choices=['ASPA','G6PD']);p.add_argument('position',type=int);p.add_argument('ref');p.add_argument('alt')
    a=p.parse_args();report=TargetFreeEffectProvider(a.bundle,a.manifest_sha256).query(a.gene,a.position,a.ref,a.alt)
    print(json.dumps(report,sort_keys=True,indent=2,allow_nan=False))


if __name__=='__main__':main()
