"""Rebuild held-residue G6PD results and verify the prediction lock using the Python standard library."""
import argparse
import csv
import gzip
import io
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median

from predict import HERE, classify, digest


def write_rows(path, rows):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_reporter(path):
    """Export the held evaluation reporter rows from the verified local provenance source."""
    expected = json.loads((HERE / "score-lock.json").read_text())["source_hashes"]["reporter"]
    if digest(path) != expected:
        raise ValueError("Reporter data hash mismatch")
    with path.open() as stream:
        for line in stream:
            row = json.loads(line)
            yield (f"{row['ref_aa']}{row['position_1based']}{row['alt_aa']}",
                   row["abundance_id"], row["abundance_score"])


def evaluate(raw=None):
    lock = json.loads((HERE / "score-lock.json").read_text())
    for name, field in [("PROTOCOL.md", "protocol_sha256"), ("model.json", "model_sha256"),
                        ("predict.py", "predictor_sha256"), ("predictions.csv.gz", "predictions_sha256")]:
        if digest(HERE / name) != lock[field]:
            raise ValueError(f"Score lock mismatch for {name}")
    model = json.loads((HERE / "model.json").read_text())
    cut = model["thresholds"]
    with gzip.open(HERE / "predictions.csv.gz","rt") as fp:
        scores = {r["substitution"]: r for r in csv.DictReader(fp)}
    assert len(scores) == lock["public_predictions"]
    with (HERE / "assessment_sources.tsv").open() as fp:
        sources = list(csv.DictReader(fp, delimiter="\t"))
    anchors = {r["source"]: float(r["kcat_s_inv"]) for r in sources if r["status"] == "anchor"}
    counts = Counter()
    direct = defaultdict(list)
    for r in sources:
        key = r["substitution"]
        status = r["status"]
        if status == "single":
            if int(key[1:-1]) in model["training_residues"]:
                status = "training_residue"
            elif key not in scores:
                status = "missing_features"
            else:
                assert scores[key]["training_residue"] == "0"
                direct[key].append((r["source"],float(r["kcat_s_inv"])/anchors[r["source"]]))
                status = "numeric_held_turnover"
        counts[status] += 1
    if raw:
        reporter = [dict(substitution=k, score_id=i, abundance=v) for k,i,v in read_reporter(raw) if k in direct and v is not None]
        write_rows(HERE / "held_reporter.csv", sorted(reporter,key=lambda r:r["substitution"]))
    with (HERE / "held_reporter.csv").open() as fp:
        reporter = list(csv.DictReader(fp))
    assert len({r["substitution"] for r in reporter}) == len(reporter)
    observed_abundance = {r["substitution"]: float(r["abundance"]) for r in reporter}
    paired = []
    for key, observations in sorted(direct.items()):
        if key not in observed_abundance:
            counts["missing_reporter"] += 1
            continue
        activity = median(v for source,v in observations)
        abundance = observed_abundance[key]
        row = scores[key]
        turnover, control, predicted_abundance = (float(row[n]) for n in ["turnover_fraction","control_turnover","reporter_abundance"])
        assert row["phenotype"] == classify(turnover,predicted_abundance,cut)
        paired.append(dict(substitution=key, residue=int(key[1:-1]), studies="|".join(s for s,v in observations),
                           observed_turnover=activity, observed_reporter=abundance, observed=classify(activity,abundance,cut),
                           predicted_turnover=turnover, predicted_reporter=predicted_abundance, predicted=row["phenotype"],
                           control_turnover=control, control=classify(control,predicted_abundance,cut)))
    def metrics(cases):
        return dict(n=len(cases),correct=sum(r["predicted"]==r["observed"] for r in cases),
                    control_correct=sum(r["control"]==r["observed"] for r in cases))
    grouped = {a+"__"+b: metrics([r for r in paired if r["observed"] == a+"__"+b])
               for a in ["low_activity","retained_activity"] for b in ["low_abundance","preserved_abundance"]}
    def mae(v,t):
        return mean(abs(r[v]-r[t]) for r in paired)
    summary=dict(coverage=dict(source_rows=len(sources), source_row_status=dict(sorted(counts.items())),
                               paired_variants=len(paired), paired_residues=len({r["residue"] for r in paired})),
                 per_observed_profile=grouped,overall=metrics(paired),
                 activity_low_binary_correct=sum((r["observed_turnover"]<cut["activity"]) == (r["predicted_turnover"]<cut["activity"]) for r in paired),
                 abundance_low_binary_correct=sum((r["observed_reporter"]<cut["abundance"]) == (r["predicted_reporter"]<cut["abundance"]) for r in paired),
                 mae=dict(turnover=mae("observed_turnover","predicted_turnover"),
                          control_turnover=mae("observed_turnover","control_turnover"),reporter=mae("observed_reporter","predicted_reporter")),
                 score_lock_sha256=digest(HERE/"score-lock.json"), assessment_sha256=digest(HERE/"assessment_sources.tsv"),
                 reporter_extract_sha256=digest(HERE/"held_reporter.csv"))
    return paired, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-reporter",type=Path,help="Re-extract paired historical abundance from the hash-verified study JSONL.")
    parser.add_argument("--write",action="store_true",help="Write result artifacts. Default verifies existing results.")
    args=parser.parse_args()
    if args.source_reporter and not args.write:
        parser.error("Re-extraction requires --write")
    rows, summary=evaluate(args.source_reporter)
    if args.write:
        write_rows(HERE/"evaluation_cases.csv",rows)
        (HERE/"results.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    else:
        if summary != json.loads((HERE/"results.json").read_text()):
            raise ValueError("Aggregate mismatch")
        buff=io.StringIO()
        writer=csv.DictWriter(buff,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        if buff.getvalue() != (HERE/"evaluation_cases.csv").read_text():
            raise ValueError("Per-variant result mismatch")
    print(json.dumps(summary, indent=2))
