"""Count historical training profiles without fitting, scoring or reading excluded variants."""
import csv
import json
from collections import Counter, defaultdict
from statistics import median

from predict import HERE, classify, digest


def inventory():
    lock = json.loads((HERE/"score-lock.json").read_text())
    model = json.loads((HERE/"model.json").read_text())
    base = HERE.parent
    turnover = base/"catalytic-effect-feasibility-v1/measurements/measurement_inventory.csv"
    reporter = base/"model-3d-understanding-v1/g6pd-paired-effects-v1/output/paired-rows.jsonl"
    for path, key in [(turnover, "turnover"), (reporter, "reporter")]:
        if digest(path) != lock["source_hashes"][key]:
            raise ValueError(f"{key} source hash mismatch")
    training = set(model["training_turnover_variants"])
    with turnover.open() as src:
        numeric = [r for r in csv.DictReader(src) if r["status"] == "numeric_available"
                   and r["endpoint"] == "kcat" and r["added_ligand"] == "not_applicable"]
    wt = {r["study"]: float(r["value"]) for r in numeric if r["background_type"] == "WT"}
    activity = defaultdict(list)
    for r in numeric:
        if r["background_type"] == "single" and r["substitutions"] in training:
            activity[r["substitutions"]].append(float(r["value"])/wt[r["study"]])
    observed = {}
    with reporter.open() as src:
        for r in map(json.loads,src):
            key = f"{r['ref_aa']}{r['position_1based']}{r['alt_aa']}"
            if key in training and r["abundance_score"] is not None:
                observed[key] = classify(median(activity[key]), r["abundance_score"], model["thresholds"])
    assert set(observed) == training
    held = json.loads((HERE/"results.json").read_text())["per_observed_profile"]
    counts = Counter(observed.values())
    return dict(training_variants=len(observed), training_labels=observed,
                training_counts={k: counts[k] for k in held},
                training_plus_primary_profile_counts={k: counts[k] + held[k]["n"] for k in held})


if __name__ == "__main__":
    print(json.dumps(inventory(), indent=2, sort_keys=True))
