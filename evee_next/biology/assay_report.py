"""Write a target-free experimental prediction from one substitution or a same-residue pair."""
import argparse
import json
from pathlib import Path

from .named_effect_provider import EFFECTS, NamedEffectProvider
from .reference_cache_provider import snapshot
from .routing import EvidenceRequest, canonical_sha, public_value
from .structure_provider import CANONICAL

BUNDLE = Path(__file__).resolve().parents[2] / 'predictor_v1/assets/effects'
PIN = '67e9757116b94e933ca8f20856d069beea004270f6d3daeb6e1aacb1c13791b1'
READOUTS = {
    'abundance': {
        'title': 'Relative reporter abundance',
        'inputs': 'Protein sequence compatibility, amino-acid chemistry, predicted stability and reference packing.',
        'limitation': 'Percentile within a matched reporter panel. The score does not quantify native protein amount.',
    },
    'growth': {
        'title': 'G6PD yeast complementation growth',
        'inputs': 'Substitution chemistry, sequence compatibility, predicted stability, reference packing, and comparison with other replacements at that site.',
        'limitation': 'Normalized yeast growth depends on several processes. The score does not measure purified enzyme activity.',
    },
}

def _case(identity, clinical, decisions):
    effects = {}
    trace = []
    for decision in decisions:
        trace.append({'provider': decision.request.provider_id, 'status': decision.status,
                      'reason': decision.reason, 'clinical_score_fields': list(decision.scoring_used)})
        for name, (field, unit, _) in EFFECTS.items():
            if decision.request.provider_id != 'effect.' + name:
                continue
            if name in effects:
                raise ValueError('Ambiguous endpoint evidence')
            record = {'status': decision.status, 'reason': decision.reason}
            if decision.status == 'available':
                if decision.units.get(field) != unit or decision.request.parameters.get('identity_status') != 'qualified':
                    raise ValueError('Unqualified endpoint identity or unit')
                for key in ('protein_id', 'residue_1based', 'ref_aa', 'alt_aa'):
                    if identity.get(key) != decision.request.identity.get(key):
                        raise ValueError('Endpoint differs from requested substitution')
                record.update(value=decision.values[field], unit=unit, qualification=decision.qualification,
                              uncertainty=decision.uncertainty, assessment=public_value(decision.metadata),
                              recipe_sha256=decision.request.recipe_sha256, provenance=public_value(decision.provenance))
            effects[name] = record
    for name in EFFECTS:
        effects.setdefault(name, {'status': 'unavailable', 'reason': 'Endpoint provider was not configured'})
    return {'identity': identity, 'clinical': clinical, 'evidence_trace': trace, 'effects': effects}

def routed_case(report):
    context = report.context
    identity = public_value(context.mapping) if context.mapping is not None else {}
    identity.update(genomic_allele=public_value(context.variant), mapping_status=context.status,
                    genomic_validated=context.genomic_validated, reference_source=public_value(context.provenance))
    clinical = {'baseline': public_value(report.baseline), 'final': public_value(report.final),
                'call': report.final.call, 'fusion_status': report.fusion_status}
    return _case(identity, clinical, report.evidence)

def protein_case(providers, manifest, protein_id, pos, ref, alt):
    refs = manifest['references']
    entry = next((r for r in refs.values() if r['protein_id'] == protein_id), None)
    if entry is None:
        raise ValueError('Versioned protein is outside the target-free packet')
    if not 1 <= pos <= len(entry['sequence']) or entry['sequence'][pos - 1] != ref:
        raise ValueError('Reference residue differs from the full versioned sequence')
    if len(alt) != 1 or alt not in CANONICAL or alt == ref:
        raise ValueError('A different canonical replacement is required')
    identity = {'organism': 'Homo sapiens', 'protein_id': protein_id, 'protein_sha256': entry['sha256'],
                'residue_1based': pos, 'ref_aa': ref, 'alt_aa': alt}
    decisions = []
    for provider in providers:
        request = EvidenceRequest(provider.provider_id, provider.recipe_sha256, identity, {'identity_status': 'qualified'})
        decisions.append(provider.obtain(request, store=None, policy=None))
    return _case(identity, {'status': 'unscored', 'reason': 'Protein-only request lacks a genomic allele'}, decisions)

def archive_case(request):
    from evee_next.report import normalize
    identity = {'genomic_allele': public_value(normalize(request)), 'mapping_status': 'paper_lookup_only',
                'requested_transcript': request.transcript_id, 'requested_protein': request.protein_id}
    return _case(identity, {'status': 'unscored', 'reason': 'Live retained clinical model was not requested'}, ())

def build_report(first, second=None):
    cases = [first] if second is None else [first, second]
    if second is not None:
        a, b = (c['identity'] for c in cases)
        same_keys = ('protein_id', 'residue_1based', 'ref_aa')
        if not all(a.get(k) and a[k] == b.get(k) for k in same_keys) or not a.get('alt_aa') or not b.get('alt_aa') or a['alt_aa'] == b['alt_aa']:
            raise ValueError('Comparison requires distinct replacements at one versioned reference residue')
    predictions = {}
    for name, (_, unit, _) in EFFECTS.items():
        records = [case['effects'][name] for case in cases]
        item = dict(READOUTS[name], unit=unit)
        if any(r['status'] != 'available' for r in records):
            item.update(status='unavailable', reasons=[r['reason'] for r in records if r['status'] != 'available'])
        else:
            for key in ('manifest_sha256', 'packet_sha256', 'protein_sha256'):
                if len({r['provenance'][key] for r in records}) != 1:
                    raise ValueError('Assay predictions use different sources or protein sequences')
            gene = records[0]['assessment']['gene']
            if name == 'growth':
                item['assay'] = 'Yeast complementation by human G6PD, normalized to the source assay controls.'
            elif gene == 'SOD1':
                item['assay'] = 'Enrichment in the high-abundance reporter selection, ranked within the SOD1 panel.'
            else:
                item['assay'] = f'Fluorescence-bin abundance of the tagged protein in human cells, ranked within the {gene} panel.'
            values = [r['value'] for r in records]
            item.update(status='predicted', values=values, uncertainty='Individual prediction uncertainty is uncalibrated.')
            if len(values) == 2:
                delta = values[0] - values[1]
                direction = 'first_higher' if delta > 0 else 'second_higher' if delta < 0 else 'tie'
                item.update(first_minus_second=delta, predicted_order=direction)
                if delta:
                    item['test'] = 'Measure both substitutions in the same assay and compare the prespecified direction.'
                    item['contradiction'] = 'A reproducible difference in the opposite direction would contradict this prediction. An unresolved difference is inconclusive.'
                else:
                    item['test'] = 'The model gives equal point estimates and makes no directional claim.'
            else:
                item['test'] = 'Compare the measured value on the source assay scale with the frozen estimate; report the signed error.'
        predictions[name] = item
    return {'schema': 'missense-assay-report/v1', 'cases': cases, 'predictions': predictions,
            'prediction_sha256': canonical_sha({'cases': cases, 'predictions': predictions})}

def render_markdown(report):
    lines = ['# Testable missense predictions', '']
    for case in report['cases']:
        i = case['identity']
        substitution = f"{i['ref_aa']}{i['residue_1based']}{i['alt_aa']}" if i.get('alt_aa') else 'Unresolved substitution'
        allele = i.get('genomic_allele')
        heading = f"{i['protein_id']} {substitution}" if i.get('protein_id') else 'Genomic request'
        if allele:
            heading += f" · {allele['assembly']} chr{allele['chromosome']}:{allele['position_1based']} {allele['ref']}>{allele['alt']} (1-based)"
        lines.append('- ' + heading)
        clinical = case['clinical'].get('final', case['clinical'])
        if clinical.get('score') is None:
            lines.append(f"  Retained clinical score unavailable. {clinical.get('reason') or clinical['status']}.")
        else:
            call = 'positive' if case['clinical']['call'] else 'negative'
            lines.append(f"  Clinical call {call}. Score {clinical['score']:.4f} at cutoff {clinical['cutoff']:.4f}. Checkpoint {clinical['checkpoint']}.")
        consumed = [row['provider'] + ' (' + ', '.join(row['clinical_score_fields']) + ')'
                    for row in case['evidence_trace'] if row['clinical_score_fields']]
        if clinical.get('score') is not None:
            if consumed:
                lines.append('  Additional evidence used in the clinical score includes ' + ', '.join(consumed) + '.')
            else:
                lines.append('  The baseline supplied the clinical score without any additional evidence fields.')
        paper = case.get('paper')
        if paper:
            if paper['status'] != 'available':
                lines.append('  Paper-wide score unavailable. ' + paper['reason'] + '.')
                continue
            call = 'positive' if paper['original_call'] else 'negative'
            lines.append(f"  The frozen paper-wide {paper['original_route']} route scored {paper['original_score']:.4f}, {call} at {paper['cutoff']:.4f}.")
            rec = paper['recovery']
            if rec['status'] == 'diagnostic_only':
                m = rec['mapping']
                lines.extend(['',
                    f"  Exploratory protein rescore {rec['score']:.4f} using {rec['route']}. The original score above remains the comparator.",
                    f"  {m['mapping_source']} · {m['transcript_id']} · {m['protein_id']} {m['ref']}{m['position']}{m['alt']} · {m['codon']}>{m['alt_codon']}."])
                if rec['alternative_transcript']:
                    lines.append(f"  The archive annotated transcript {m['archived_transcript']} as {m['annotation']}. Verify expression of the alternative isoform in the disease-relevant tissue before interpreting the protein score.")
            elif rec['status'] == 'unavailable':
                lines.append('  Protein rescore unavailable. ' + rec['reason'] + '.')
    for name, item in report['predictions'].items():
        lines.extend(['', '## ' + item['title'], ''])
        if item['status'] != 'predicted':
            lines.append('Unavailable. ' + '. '.join(dict.fromkeys(item['reasons'])).rstrip('.') + '.')
            continue
        for case, value in zip(report['cases'], item['values']):
            i = case['identity']
            display = f"reporter abundance percentile {100*value:.1f}" if name == 'abundance' else f'{value:.3f} normalized growth'
            lines.append(f"- {i['ref_aa']}{i['residue_1based']}{i['alt_aa']} — {display}.")
        if 'predicted_order' in item and item['predicted_order'] != 'tie':
            high = report['cases'][0 if item['predicted_order'] == 'first_higher' else 1]['identity']
            lines.extend(['', f"The model predicts higher {item['title'].replace('Relative', 'relative')} for {high['ref_aa']}{high['residue_1based']}{high['alt_aa']}."])
        model_use = report['cases'][0]['effects'][name]['assessment']
        if model_use['assessment'] == 'held_out_protein':
            assessment = ('The model was fitted on ' + ', '.join(model_use['fitting_proteins']) +
                          '. The assessment held out all substitutions in ' + model_use['gene'] + '.')
        else:
            start = model_use['held_region_0based'] * 103 + 1
            assessment = f'The G6PD assessment held out the region containing residues {start}–{start+102}.'
        endpoint_used = any(name_record['clinical_score_fields'] for case in report['cases']
                            for name_record in case['evidence_trace'] if name_record['provider'] == 'effect.' + name)
        contribution = 'The clinical scoring route consumed this endpoint.' if endpoint_used else 'This endpoint did not contribute to a clinical score.'
        lines.extend(['', item['assay'], item['test'], item.get('contradiction', ''), item['uncertainty'], item['limitation'],
                      '', assessment, contribution,
                      '', 'Model inputs include ' + item['inputs'][:1].lower() + item['inputs'][1:]])
        provenance = report['cases'][0]['effects'][name]['provenance']
        lines.extend(['', 'Prediction packet SHA-256 `' + provenance['packet_sha256'] + '`.'])
    return '\n'.join(lines) + '\n'

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('protein', help='Versioned Ensembl protein ID')
    parser.add_argument('position', type=int)
    parser.add_argument('ref')
    parser.add_argument('alt')
    parser.add_argument('--compare-alt')
    parser.add_argument('--bundle', type=Path, default=BUNDLE)
    parser.add_argument('--manifest-sha256', default=PIN)
    parser.add_argument('--format', choices=('json', 'markdown'), default='json')
    args = parser.parse_args(argv)
    try:
        manifest = json.loads(snapshot(args.bundle/'manifest.json', args.manifest_sha256, 64*1024))
        providers = [NamedEffectProvider(args.bundle, args.manifest_sha256, name) for name in EFFECTS]
        def case(alt):
            return protein_case(providers, manifest, args.protein, args.position, args.ref, alt)
        report = build_report(case(args.alt), case(args.compare_alt) if args.compare_alt else None)
    except (KeyError, OSError, ValueError) as exc:
        parser.exit(2, f'{exc}\n')
    print(render_markdown(report) if args.format == 'markdown' else json.dumps(report, indent=2, allow_nan=False), end='\n')

if __name__ == '__main__':
    main()
