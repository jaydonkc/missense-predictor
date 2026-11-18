"""Query exact held-out human assay predictions from target-free development packets."""
import json
import math
from pathlib import Path
from types import MappingProxyType

from .reference_cache_provider import snapshot
from .routing import EvidenceDecision, canonical_bytes, canonical_sha, require_sha
from .structure_provider import CANONICAL, _digest, _request

# Units and the permitted interpretation are fixed by the original endpoint.
EFFECTS = {
    'abundance': ('relative_reporter_abundance', 'within_assay_percentile_fraction_0_to_1', 'held_out_protein'),
    'growth': ('yeast_complementation_growth', 'normalized_yeast_complementation_growth', 'held_out_G6PD_region'),
}
PROTEINS = {'ASPA': 'ENSP00000263080.2', 'G6PD': 'ENSP00000377192.3', 'SOD1': 'ENSP00000270142.7'}


class NamedEffectProvider:
    def __init__(self, bundle, manifest_sha256, endpoint):
        self.endpoint = endpoint
        self.field, self.unit, self.assessment = EFFECTS[endpoint]
        self.provider_id = 'effect.' + endpoint
        self.pin = manifest_sha256
        folder = Path(bundle)
        manifest = json.loads(snapshot(folder/'manifest.json', self.pin, 64*1024))
        if (set(manifest) != {'schema','references','endpoints'} or
                manifest['schema'] != 'human-missense-named-effects/v1' or
                set(manifest['endpoints']) != set(EFFECTS)):
            raise ValueError('Named effect packet schema mismatch')
        spec = manifest['endpoints'][endpoint]
        if (set(spec) != {'field','unit','assessment','arm','packet_sha256','source_predictions_sha256',
                          'source_verification_sha256','rows'} or
                (spec['field'], spec['unit'], spec['assessment']) != EFFECTS[endpoint] or
                spec['arm'] != {'growth':'nonlinear','abundance':'candidate'}[endpoint]):
            raise ValueError('Unqualified effect source or units')
        for k in ('packet_sha256','source_predictions_sha256','source_verification_sha256'):
            require_sha(spec[k])
        self.spec = MappingProxyType(spec)
        refs = manifest['references']
        if set(refs) != set(PROTEINS):
            raise ValueError('Expected three exact human protein references')
        self.refs = {}
        for gene, r in refs.items():
            if (set(r) != {'protein_id','sequence','sha256'} or r['protein_id'] != PROTEINS[gene] or
                    not r['sequence'] or not set(r['sequence']) <= CANONICAL or
                    _digest(r['sequence'].encode()) != r['sha256']):
                raise ValueError('Inconsistent human protein reference')
            self.refs[r['protein_id']] = (gene, r['sha256'])
        rows = json.loads(snapshot(folder/(endpoint+'.json'), spec['packet_sha256'], 8*1024**2))
        if not isinstance(rows, list) or not rows or type(spec['rows']) is not int or len(rows) != spec['rows']:
            raise ValueError('Effect roster size differs')
        self.lookup = {}
        for gene, pos, ref, alt, value, fold in rows:
            seq = refs[gene]['sequence']
            if (type(pos) is not int or not 1 <= pos <= len(seq) or seq[pos-1] != ref or
                    alt == ref or alt not in CANONICAL or type(value) not in (float,int) or not math.isfinite(value)):
                raise ValueError('Invalid effect roster substitution or value')
            if endpoint == 'abundance' and (fold is not None or not 0 <= value <= 1):
                raise ValueError('Abundance percentile or holdout differs')
            if endpoint == 'growth' and (gene != 'G6PD' or type(fold) is not int or fold != (pos-1)//103):
                raise ValueError('Growth source protein or held region differs')
            key = (gene,pos,ref,alt)
            if key in self.lookup:
                raise ValueError('Duplicate effect roster substitution')
            self.lookup[key] = (value,fold)
        self.recipe_sha256 = canonical_sha({'provider':self.provider_id,'manifest':self.pin,
                                            'code_sha256':_digest(Path(__file__).read_bytes()),
                                            'identity_code_sha256':_digest(Path(__file__).with_name('structure_provider.py').read_bytes())})

    def plan(self, context):
        return (_request(self.provider_id, self.recipe_sha256, context),)

    def obtain(self, request, store, policy):
        if request.provider_id != self.provider_id or request.recipe_sha256 != self.recipe_sha256:
            return EvidenceDecision(request, True, 'failed', 'Named effect source/recipe mismatch')
        def produce():
            i = request.identity
            protein = self.refs.get(i.get('protein_id'))
            if request.parameters.get('identity_status') != 'qualified' or i.get('organism') != 'Homo sapiens':
                return b'', EvidenceDecision(request, True, 'unavailable', 'A unique, qualified human protein mapping is required')
            if protein is None or i.get('protein_sha256') != protein[1]:
                return b'', EvidenceDecision(request, True, 'unavailable', 'No matching versioned full protein sequence in this effect packet')
            gene = protein[0]
            if self.endpoint == 'growth' and gene != 'G6PD':
                return b'', EvidenceDecision(request, True, 'unavailable', 'Growth model is qualified only for human G6PD')
            value = self.lookup.get((gene,i['residue_1based'],i['ref_aa'],i['alt_aa']))
            if value is None:
                return b'', EvidenceDecision(request, True, 'unavailable', 'Exact amino acid substitution absent from the effect packet')
            score, fold = value
            metadata = {'gene':gene, 'assessment':self.assessment,'source_arm':self.spec['arm'],
                        'prediction_type':'assay_effect_prediction','clinical_use':'diagnostic_only'}
            if self.endpoint == 'growth': metadata['held_region_0based'] = fold
            else: metadata['fitting_proteins'] = {'ASPA':['G6PD'], 'G6PD':['ASPA'], 'SOD1':['ASPA','G6PD']}[gene]
            provenance = {'manifest_sha256':self.pin,'protein_sha256':protein[1],
                          **{k:self.spec[k] for k in ('packet_sha256','source_predictions_sha256','source_verification_sha256')}}
            d = EvidenceDecision(request,True,'available','Exact cached assay-model prediction',
                {self.field:score},{self.field:self.unit},provenance,metadata,qualification='development_diagnostic',
                uncertainty='Individual prediction uncertainty is uncalibrated. Evaluation reused exposed research panels.',
                execution='local_named_effect_lookup')
            return canonical_bytes(d), d
        if store is None: return produce()[1]
        return store.get_or_compute(request, produce)
