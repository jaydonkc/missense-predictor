"""Render a target-free report for a frozen human G6PD prediction."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def verified(path, expected):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Frozen artifact failed its integrity check: {path.name}")
    return path


def report(substitution):
    lock = json.loads((HERE / "score-lock.json").read_text())
    state = json.loads(verified(HERE / "model.json", lock["model_sha256"]).read_bytes())
    path = verified(HERE / "predictions.csv.gz", lock["predictions_sha256"])
    with gzip.open(path, "rt") as fp:
        row = next((r for r in csv.DictReader(fp) if r["substitution"] == substitution), None)
    if row is None:
        raise ValueError(f"No qualified prediction for {substitution} on {state['protein_id']}")
    return {
        "gene": "G6PD", "protein": state["protein_id"], "substitution": substitution,
        "turnover_fraction": float(row["turnover_fraction"]),
        "reporter_abundance": float(row["reporter_abundance"]),
        "thresholds": state["thresholds"],
        "predicted_assay_profile": row["phenotype"],
        "training_residue": row["training_residue"] == "1",
        "checkpoint": lock["model_sha256"],
        "role": "exploratory diagnostic; excluded from clinical scoring",
        "uncertainty": "Individual prediction uncertainty is uncalibrated.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("substitution", help="Single amino-acid change on the G6PD reference, e.g. T279A")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    try:
        result = report(args.substitution)
    except ValueError as error:
        parser.error(str(error))
    if args.json:
        print(json.dumps(result, indent=2, allow_nan=False))
        return
    thresholds = result["thresholds"]
    print(f"## {result['gene']} {result['substitution']} on {result['protein']}")
    print(f"\nPredicted catalytic turnover is {result['turnover_fraction']:.3f} times matched wild type. "
          f"Predicted human-cell reporter abundance is {result['reporter_abundance']:.3f} in normalized assay units.")
    activity, abundance = result["predicted_assay_profile"].split("__")
    print(f"\nAssay hypothesis: {activity.replace('_', ' ')} with {abundance.replace('_', ' ')}.")
    print(f"\nTest the two predictions independently. Compare purified mutant and wild-type turnover in the same "
          f"assay and measure normalized abundance with the matched human-cell reporter. The frozen boundaries "
          f"are {thresholds['activity']} for turnover and {thresholds['abundance']} for abundance. "
          "A reproducible result across either boundary from its prediction contradicts that part of the hypothesis.")
    print("\nThis exploratory model has no effect on the clinical call. " + result["uncertainty"])
    if result["training_residue"]:
        print("This position was used during endpoint fitting.")

if __name__ == "__main__":
    main()
