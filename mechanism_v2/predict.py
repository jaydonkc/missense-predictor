"""Predict normalized catalytic turnover and reporter abundance for human G6PD."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def linear(head, values):
    if len(values) != len(head["coef"]) or any(not math.isfinite(v) for v in values):
        raise ValueError("Supply a finite value for each required molecular feature")
    return head["intercept"] + sum(c * x for c, x in zip(head["coef"], values))


def classify(activity, abundance, thresholds):
    first = "low_activity" if activity < thresholds["activity"] else "retained_activity"
    second = "low_abundance" if abundance < thresholds["abundance"] else "preserved_abundance"
    return first + "__" + second


def predict(state, compatibility, sequence_likelihood, stability, abundance_percentile):
    activity = max(0., linear(state["activity"], [compatibility, sequence_likelihood, stability]))
    abundance = linear(state["abundance"], [abundance_percentile])
    control = max(0., linear(state["activity_control"], [abundance_percentile]))
    return dict(turnover_fraction=activity, reporter_abundance=abundance,
                phenotype=classify(activity, abundance, state["thresholds"]),
                control_turnover=control)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--variant", help="Human G6PD substitution in one-letter format")
    mode.add_argument("--features", nargs=4, type=float, metavar=("SAPROT", "ESMC", "DDG", "ABUNDANCE_PERCENTILE_0_TO_1"))
    args = parser.parse_args()
    lock = json.loads((HERE/"score-lock.json").read_text())
    state_path = HERE/"model.json"
    if digest(state_path) != lock["model_sha256"]:
        raise ValueError("Model hash mismatch")
    state = json.loads(state_path.read_text())
    if args.features is not None:
        if not 0 <= args.features[3] <= 1:
            parser.error("Abundance percentile fraction must be between zero and one")
        answer = predict(state, *args.features)
    else:
        path = HERE/"predictions.csv.gz"
        if digest(path) != lock["predictions_sha256"]:
            raise ValueError("Prediction hash mismatch")
        with gzip.open(path, "rt") as stream:
            answer = next((r for r in csv.DictReader(stream) if r["substitution"] == args.variant), None)
        if answer is None:
            parser.error("No qualified prediction for this substitution on " + state["protein_id"])
    print(json.dumps(dict(protein=state["protein_id"], **answer), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
