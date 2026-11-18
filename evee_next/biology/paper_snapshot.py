"""Read exact genomic-allele predictions from two locked, outcome-free paper exports."""
import csv
import gzip
import io
import json

from evee_next.report import normalize
from .reference_cache_provider import snapshot

PAPER_PIN = '97c1f33fb781bd9f36392c67509b0368db6ce2723a5738db2538c4782ee25d89'
RECOVERY_PIN = 'ea528b1817ab3be0f0a59435eaa9fdd9f1e6f477db06acfdff08413d4f4d3253'

def _selected(path, pin, ids):
    data = snapshot(path, pin, 8*1024**2)
    rows = {}
    with gzip.open(io.BytesIO(data), 'rt', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            key = row['variant_id']
            if key in ids:
                if key in rows:
                    raise ValueError('Duplicate paper allele')
                rows[key] = row
    return rows

def lookups(study, requests):
    ids = {normalize(request).variant_id for request in requests}
    recovered_dir = study/'paper-protein-recovery-v1'
    paper = json.loads(snapshot(study/'paper-cutoff-transfer-v1/score-lock.json', PAPER_PIN))
    recovery = json.loads(snapshot(recovered_dir/'score-lock.json', RECOVERY_PIN))
    original = _selected(study/'paper-full-missense-v1/paper-scores.tsv.gz', paper['inputs']['scores']['sha256'], ids)
    replacements = _selected(recovered_dir/'recovered-scores.tsv.gz', recovery['scores_sha256'], ids)
    mappings = _selected(recovered_dir/'mappings.tsv.gz', recovery['files']['mappings.tsv.gz'], ids)
    ledger = _selected(recovered_dir/'fallback-ledger.tsv.gz', recovery['files']['fallback-ledger.tsv.gz'], ids)
    if replacements.keys() != mappings.keys():
        raise ValueError('Recovery prediction and molecular mapping differ')
    result = {}
    for key in ids:
        record = {'allele': key, 'scope': 'exact GRCh38 genomic allele from the archived paper cohort',
                  'paper_lock_sha256': PAPER_PIN, 'recovery_lock_sha256': RECOVERY_PIN}
        row = original.get(key)
        if row is None:
            record.update(status='unavailable', reason='No exact allele in frozen paper score export')
            result[key] = record
            continue
        score, cutoff = float(row['score']), paper['cutoff']
        record.update(status='available', model=paper['model'], original_route=row['route'], original_score=score,
                      cutoff=cutoff, original_call=score >= cutoff)
        if key in replacements:
            recovered, mapping = replacements[key], mappings[key]
            if (row['route'] != 'DNA only' or float(recovered['previous_score']) != score or
                    mapping['mapping_source'] != recovered['mapping_source']):
                raise ValueError('Recovery does not reconcile with the original DNA score')
            source = mapping['mapping_source']
            record['recovery'] = dict(status='diagnostic_only', mapping_source=source, mapping=mapping,
                alternative_transcript=source != 'paper transcript', score=float(recovered['score']),
                route=recovered['route'], above_original_cutoff=float(recovered['score']) >= cutoff)
        elif key in ledger:
            record['recovery'] = {'status': 'unavailable', 'reason': ledger[key]['disposition']}
        elif row['route'] == 'DNA only':
            raise ValueError('DNA-only allele missing its recovery disposition')
        else:
            record['recovery'] = {'status': 'not_applicable', 'reason': 'Recovery targeted the DNA-only route'}
        result[key] = record
    return result
