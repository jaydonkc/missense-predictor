"""Replay the frozen public comparison. Run from any directory with Python and NumPy."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))


def read_table(name, manifest):
    path = HERE / name
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest['files'][name]:
        raise ValueError(f'Checksum mismatch for {name}')
    with gzip.open(path, 'rt', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def development_cutoff(rows, model):
    """Minimize FP + FN, keeping ties together and preferring proximity to 0.5."""
    counts = {float(r['score']): (int(r['pathogenic']), int(r['benign']))
              for r in rows if r['model'] == model}
    misses = sum(pos for pos, neg in counts.values())
    false_alarms = 0
    assert misses == 7481 and sum(neg for pos, neg in counts.values()) == 17867
    best = None
    candidates = set(counts) | {0.5, math.nextafter(max(counts), math.inf)}
    for cutoff in sorted(candidates, reverse=True):
        pos, neg = counts.get(cutoff, (0, 0))
        misses -= pos
        false_alarms += neg
        key = (misses + false_alarms, abs(cutoff - .5), -cutoff)
        if best is None or key < best:
            best = key
    return -best[2]


def paired_summary(ours, evee, labels, genes):
    our_errors, other_errors = ours != labels, evee != labels
    _, groups = np.unique(genes, return_inverse=True)
    sizes = np.bincount(groups)
    net = np.bincount(groups, weights=other_errors.astype(int) - our_errors.astype(int))
    rng = np.random.default_rng(1729)
    draws = []
    for _ in range(50):
        sampled = rng.integers(len(sizes), size=(100, len(sizes)))
        draws.extend(1000 * net[sampled].sum(axis=1) / sizes[sampled].sum(axis=1))
    return dict(repaired=int(sum(other_errors & ~our_errors)),
                introduced=int(sum(our_errors & ~other_errors)),
                fewer_errors_per_1000_ci95=np.quantile(draws, [.025, .975]).round(3).tolist())


def replay(prepared):
    manifest = json.loads((HERE / 'manifest.json').read_text())
    rows = read_table('scored_variants.tsv.gz', manifest)
    dev = read_table('development_counts.tsv.gz', manifest)
    ids = [r['variant_id'] for r in rows]
    assert len(ids) == len(set(ids)) == 10000
    evidence = read_table('evidence_availability.tsv.gz', manifest)
    assert len(evidence) == len(ids) and {r['variant_id'] for r in evidence} == set(ids)
    assert all(r['esmc_available'] == 'True' for r in evidence)
    assert sum(r['reference_structure'] == 'observed' for r in evidence) == 7400
    assert sum(r['burial_available'] == '1' for r in evidence) == 3943
    folds = np.array([[float(r[f'fold{i}']) for i in range(3)] for r in rows])
    if prepared is not None:
        from predictor_v1.score import Predictor
        model = Predictor()
        with np.load(prepared, allow_pickle=False) as features:
            assert features['ids'].tolist() == ids
            derived = np.column_stack([
                model.full_context_scores(i, features['adjusted_parent'][:, i], features['x'])
                for i in range(3)])
        np.testing.assert_allclose(derived, folds, rtol=0, atol=1e-12)
        folds = derived
    ours = []
    for row, scores in zip(rows, folds):
        eligible = [int(f) for f in row['eligible_folds'].split(',')]
        assert eligible in ([0], [1], [2], [0, 1, 2])
        ours.append(math.fsum(scores[f] for f in eligible) / len(eligible))
    np.testing.assert_allclose(ours, [float(r['score']) for r in rows], rtol=0, atol=1e-12)
    ours = np.asarray(ours)
    evee = np.array([float(r['evee']) for r in rows])
    y = np.array([int(r['pathogenic']) for r in rows])
    genes = [r['gene_id'] for r in rows]
    assert set(y) == {0, 1} and sum(y) == 2861 and len(set(genes)) == 3813
    a, b = (development_cutoff(dev, name) for name in ('predictor', 'EVEE'))
    result = []
    for model, score, cutoff, rule in [
        ('Our predictor', ours, a, 'Same development set'), ('EVEE', evee, b, 'Same development set'),
        ('Our predictor', ours, .5, 'Paper cutoff'), ('EVEE', evee, .5, 'Paper cutoff'),
        ('EVEE', evee, .3738, 'Prespecified comparator')]:
        pred = score >= cutoff
        fn, fp = int(sum((y == 1) & ~pred)), int(sum((y == 0) & pred))
        result.append(dict(rule=rule, model=model, cutoff=cutoff, misses=fn, false_alarms=fp,
                           errors=fn + fp, accuracy=round(1 - (fn + fp)/len(y), 4)))
    return dict(variants=len(rows), operating_points=result,
                same_development=paired_summary(ours >= a, evee >= b, y, genes),
                prespecified=paired_summary(ours >= a, evee >= .3738, y, genes))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, help='Optional original DGX intermediates.npz in the locked row order')
    args = parser.parse_args()
    actual = replay(args.prepared)
    expected = json.loads((HERE / 'results.json').read_text())
    if actual != expected:
        raise SystemExit(f'Benchmark differs from expected results: {actual}')
    print(json.dumps(actual, indent=2))
