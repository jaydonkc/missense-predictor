"""Replay the locked categorical ASPA transfer check using only the standard library."""
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CLASSES = ('selective', 'dual', 'neither')

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def evaluate():
    lock = json.loads((HERE/'score-lock.json').read_text())
    packet = HERE/'predictions.tsv.gz'
    assert digest(packet)==lock['predictions_sha256'] and digest(HERE/'PROTOCOL.md')==lock['protocol_sha256']
    with gzip.open(packet,'rt') as fp:
        predictions={r['substitution']:r for r in csv.DictReader(fp,delimiter='\t')}
    assert len(predictions)==lock['rows']
    output={}
    with (HERE/'measurements.tsv').open() as fp:
        grouped=defaultdict(list)
        for m in csv.DictReader(fp,delimiter='\t'):
            grouped[m['study']].append(m)
    for study, measurements in grouped.items():
        by_class=Counter()
        excluded=[]
        rows=[]
        confusion={arm:{c:Counter() for c in CLASSES} for arm in ['category','abundance_only']}
        assert len({r['substitution'] for r in measurements})==len(measurements)
        for m in measurements:
            key=m['substitution']
            if 'censored' in (m['activity'],m['amount']):
                excluded.append([key,'source_censored']);continue
            if key not in predictions or predictions[key]['category']=='unavailable':
                excluded.append([key,'score_unavailable']);continue
            obs={('low','preserved'):'selective',('low','low'):'dual',('preserved','preserved'):'neither'}[(m['activity'],m['amount'])]
            p=predictions[key]
            by_class[obs]+=1
            for arm in confusion:
                confusion[arm][obs][p[arm]]+=1
            rows.append(dict(substitution=key,observed=obs,prediction=p['category'],abundance=float(p['abundance_percentile']),distance_A=float(p['contact_distance_A'])))
        scores={}
        for arm,conf in confusion.items():
            scores[arm]=dict(correct=sum(conf[c][c] for c in by_class),
                             correct_by_class={c:conf[c][c] for c in by_class},
                             observed_class_macro_recall=sum(conf[c][c]/by_class[c] for c in by_class)/len(by_class))
        output[study]=dict(all_rows=len(measurements),scored_rows=len(rows),observed=by_class,excluded=excluded,
                           confusion=confusion,scores=scores,rows=rows)
    return dict(predictions_sha256=digest(packet), measurements_sha256=digest(HERE/'measurements.tsv'),studies=output)

if __name__=='__main__':
    print(json.dumps(evaluate(),indent=2))
