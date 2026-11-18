"""Diagnostic chemistry of a human missense allele at a UniProt disulfide endpoint."""
import json
from pathlib import Path

from .annotation_provider import FunctionalAnnotationProvider
from .execution import BudgetExceeded
from .routing import EvidenceDecision, canonical_sha
from .structure_provider import CANONICAL, _digest


def exact_interval(feature, length):
    """Read canonical exact endpoints only. UniProt isoform locations carry a sequence."""
    loc = feature.get('location', {})
    if set(loc) != {'start', 'end'}:
        return None
    start, end = loc['start'], loc['end']
    if any(not isinstance(x, dict) or x.get('modifier') != 'EXACT'
           or type(x.get('value')) is not int for x in (start, end)):
        return None
    a, b = start['value'], end['value']
    return (a, b) if 1 <= a <= b <= length else None


def endpoint_annotations(features, sequence, position):
    """Select simple same-chain reference bonds. A missing annotation is inconclusive."""
    if not isinstance(features, list) or any(not isinstance(f, dict) for f in features):
        raise ValueError('Native UniProt feature collection is malformed')
    chains = [exact_interval(f, len(sequence)) for f in features if f.get('type') == 'Chain']
    matches = []
    for f in features:
        if f.get('type') != 'Disulfide bond':
            continue
        pair = exact_interval(f, len(sequence))
        if pair is None or position not in pair:
            continue
        a, b = pair
        # Sites with descriptive qualifiers (interchain, alternatives, redox
        # states etc.) require separate interpretation; never reduce them to zeros.
        if (a == b or f.get('description', '').strip().lower() not in {'', 'intrachain'}
                or sequence[a-1] != 'C' or sequence[b-1] != 'C'
                or not any(c and c[0] <= a < b <= c[1] for c in chains)):
            continue
        evidence = f.get('evidences', [])
        if not isinstance(evidence, list) or any(not isinstance(e, dict) for e in evidence):
            raise ValueError('Malformed UniProt feature evidence')
        citations = [{k: e[k] for k in ('evidenceCode', 'source', 'id') if k in e} for e in evidence]
        codes = {e.get('evidenceCode') for e in citations}
        if 'ECO:0000269' in codes:
            basis = 'experimental_evidence_code_on_reference_annotation'
        elif 'ECO:0000255' in codes:
            basis = 'sequence_model_inference_on_reference_annotation'
        else:
            basis = 'other_or_unspecified_reference_evidence'
        matches.append({'endpoints_1based': [a, b], 'partner_1based': b if a == position else a,
                        'feature_id': f.get('featureId'), 'evidence_basis': basis, 'evidences': citations})
    return matches


class DisulfideLossProvider(FunctionalAnnotationProvider):
    provider_id = 'annotation.disulfide'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.recipe_sha256 = canonical_sha({'annotation_recipe_sha256': self.recipe_sha256,
                                           'disulfide_code_sha256': _digest(Path(__file__).read_bytes())})

    def _failure(self, request, status, reason):
        return EvidenceDecision(request, status != 'not_applicable', status, reason,
                                metadata={'endpoint_state': 'unresolved'})

    def _derive(self, request, store, policy):
        identity = request.identity
        state = request.parameters.get('identity_status')
        if state != 'qualified':
            return b'', self._failure(request, 'ambiguous' if state == 'ambiguous' else 'unavailable',
                                      'Unique validated missense and full protein reference required')
        if identity.get('organism') != 'Homo sapiens':
            return b'', self._failure(request, 'not_applicable', 'Route is qualified for human proteins')
        if identity.get('ref_aa') != 'C' or identity.get('alt_aa') not in CANONICAL - {'C'}:
            return b'', self._failure(request, 'not_applicable', 'Substitution does not replace a cysteine')
        records = [r for r in self.records if r.get('protein_id') == identity['protein_id']
                   and r.get('protein_sha256') == identity['protein_sha256']]
        if len(records) != 1:
            return b'', self._failure(request, 'ambiguous' if records else 'unavailable',
                                      'Expected one exact full-protein source record')
        record = records[0]
        try:
            raw, disposition, execution = self._payload(record, request, store, policy)
            data = json.loads(raw)
            if (_digest(raw) != record['sha256'] or data['primaryAccession'] != record['accession']
                    or data['entryAudit']['entryVersion'] != record['entry_version']
                    or data['organism']['taxonId'] != 9606):
                raise ValueError('Native human UniProt source identity or digest differs')
            sequence, pos = data['sequence']['value'], identity['residue_1based']
            if (_digest(sequence.encode()) != identity['protein_sha256'] or
                    not 1 <= pos <= len(sequence) or sequence[pos-1] != 'C'):
                return raw, self._failure(request, 'unavailable', 'Full protein or focal residue mismatch')
            matches = endpoint_annotations(data['features'], sequence, pos)
            provenance = {'source_id': self.manifest['source_id'], 'source_version': self.manifest['source_version'],
                          'license': self.manifest['license'], 'accession': record['accession'],
                          'entry_version': str(record['entry_version']), 'artifact_sha256': record['sha256'],
                          'manifest_sha256': self.manifest_sha256, 'protein_sha256': identity['protein_sha256']}
            partners = {m['partner_1based'] for m in matches}
            if len(partners) != 1:
                return raw, EvidenceDecision(request, True, 'ambiguous' if partners else 'unavailable',
                    'Multiple qualified reference partners' if partners else 'No qualified same-chain endpoint annotation',
                    provenance=provenance, metadata={'reference_annotations': matches},
                    cache_disposition=disposition, execution=execution)
            return raw, EvidenceDecision(request, True, 'available',
                'The replacement lacks the cysteine thiol needed for this annotated reference pair',
                {'removes_annotated_endpoint_cysteine': 1}, {'removes_annotated_endpoint_cysteine': 'binary_indicator'},
                provenance, {'reference_annotations': matches, 'partner_1based': next(iter(partners)),
                             'effect_scope': 'substitution_chemistry_given_reference_pair'},
                uncertainty='Reference pairing has the source evidence shown. Mutant folding and function are unmeasured.',
                cache_disposition=disposition, execution=execution)
        except PermissionError as exc:
            return b'', self._failure(request, 'pending_authorization', str(exc))
        except (FileNotFoundError, BudgetExceeded) as exc:
            return b'', self._failure(request, 'unavailable', str(exc))
        except (KeyError, ValueError, TypeError, OSError) as exc:
            return b'', self._failure(request, 'failed', f'Disulfide annotation extraction failed: {exc}')
