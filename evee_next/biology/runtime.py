"""Read-only intermediate-effect queries using the released cached inputs.

Protein requests do not invent a genomic allele or produce a clinical score.
For genomic requests use BiologyPipeline with the existing ReferenceResolver
and standalone DNA/protein analyzer. No training targets are read by this module.
"""
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from .chemistry import Chemistry, ProteinVariant
from .cache import load_inputs, digest, checked
from .effects import designs, predict_ridge

STUDY=Path('research/evee-next-2026-11-05/biology-pipeline-v1')
PROTEINS={'ASPA':'ENSP00000263080.2','G6PD':'ENSP00000377192.3'}


def protein_report(variant,sequence,chemistry,provider):
    variant.validate(sequence)
    effect=provider(variant) if provider else None
    return dict(schema='standalone-biology-protein-report/v1',variant=asdict(variant),
                chemistry=chemistry.calculate(variant.ref_aa,variant.alt_aa),
                molecular_effect=effect,effect_status='cached_development_effect' if effect else 'outside_cached_effect_roster',
                clinical_score=None,clinical_score_reason='protein_assay_request_has_no_resolved_genomic_DNA_inputs',
                clinical_predictor_changed=False,mechanism_status='unresolved')


class CachedEffectProvider:
    def __init__(self,root: Path):
        self.root=root;self.study=root/STUDY
        release=json.loads((self.study/'RELEASE.json').read_bytes())
        for path,sha in release['source_hashes'].items():checked(root/path,sha)
        output=self.study/'output-v1';manifest=json.loads((output/'manifest.json').read_bytes())
        # Only model/result files are needed; do not open saved target/prediction arrays.
        for key in ('heads.json','result.json'):checked(output/key,manifest[key])
        self.result=json.loads((output/'result.json').read_bytes())
        self.heads=json.loads((output/'heads.json').read_bytes())
        self.data=load_inputs(root,release,targets=False)
        self.designs,_,_=designs(self.data['chemistry'],self.data['context'],self.data['protein'],self.data['rows'])
        self.lookup={(r['gene'],r['position0']+1,r['ref_aa'],r['alt_aa']):i for i,r in enumerate(self.data['rows'])}

    def __call__(self,v):
        if PROTEINS.get(v.gene)!=v.protein_id:raise ValueError('Unsupported exact protein reference')
        if self.data['references'][v.gene]['sha256']!=v.reference_sha256:raise ValueError('Source reference differs')
        key=(v.gene,v.position_1based,v.ref_aa,v.alt_aa);i=self.lookup.get(key)
        if i is None:return None
        values={}
        for arm in ('C','R','K','P','J'):
            head,=[h for h in self.heads if h['assessment_gene']==v.gene and h['arm']==arm]
            values[arm]=float(predict_ridge(head['state'],self.designs[arm][i:i+1])[0])
        ctx=self.data['context'][i];available=bool(ctx[2])
        return dict(endpoint='reporter_abundance',units='fitting_assay_native_fluorescence_bin_scale',
                    fitting_gene='G6PD' if v.gene=='ASPA' else 'ASPA',predictions=values,
                    context=dict(reference_pLDDT=float(ctx[0]*100),reference_RSA=float(ctx[1]) if available else None,
                                 RSA_available=available,source='pinned_reference_context_cache'),
                    qualification='development_lead' if self.result['biological_development_lead'] else 'failed_advancement_do_not_use_for_clinical_scoring',
                    uncertainty='No calibrated individual prediction interval; see held-protein errors and paired development intervals',
                    source=str((self.study/'output-v1/result.json').resolve()),
                    result_sha256=digest(self.study/'output-v1/result.json'),
                    is_pathogenicity_probability=False)

    def query(self,gene,position,ref,alt):
        if gene not in PROTEINS:raise ValueError('No cached effect input for that protein')
        reference=self.data['references'][gene]
        v=ProteinVariant(gene,PROTEINS[gene],position,ref,alt,reference['sha256'])
        return protein_report(v,reference['sequence'],Chemistry.from_project(self.root),self)
