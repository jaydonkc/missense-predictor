"""Freeze V2 on historical turnover variants and their residues' abundance measurements."""
import csv
import gzip
import io
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from predict import HERE, digest, predict

STUDY = HERE.parent
ABUNDANCE = STUDY/"abundance-effect-predictor-v1/run-v1"
SOURCES = {
    "turnover": STUDY/"catalytic-effect-feasibility-v1/measurements/measurement_inventory.csv",
    "reporter": STUDY/"model-3d-understanding-v1/g6pd-paired-effects-v1/output/paired-rows.jsonl",
    "features": ABUNDANCE/"feature-cache.npz",
    "endpoint": STUDY/"named-effect-routing-v1/bundle-v1/abundance.json",
    "reference": STUDY/"named-effect-routing-v1/bundle-v1/manifest.json",
}


def fit(x, y, weights):
    x, y, w = np.asarray(x, float), np.asarray(y, float), np.asarray(weights, float)
    w = w / w.mean()
    mx, my = np.average(x, weights=w, axis=0), np.average(y, weights=w)
    scale = np.sqrt(np.average((x-mx)**2, weights=w, axis=0))
    scale[scale == 0] = 1
    z = (x - mx) / scale
    coef = np.linalg.solve(z.T @ (w[:, None]*z) + 10*np.eye(x.shape[1]), z.T @ (w*(y-my))) / scale
    return dict(intercept=float(my-mx@coef), coef=coef.tolist())


def write_new(path, data):
    with path.open("x") as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write("\n")


def train():
    ref = json.loads(SOURCES["reference"].read_text())
    protein = ref["references"]["G6PD"]
    manifest = json.loads((ABUNDANCE/"manifest.json").read_text())
    assert digest(SOURCES["features"]) == manifest["files"]["feature-cache.npz"]
    assert digest(SOURCES["endpoint"]) == ref["endpoints"]["abundance"]["packet_sha256"]
    features = {}
    with np.load(SOURCES["features"], allow_pickle=False) as z:
        for g,i,r,a,x in zip(*(z[n] for n in ["genes","positions0","ref_aa","alt_aa","X"]), strict=True):
            if g == "G6PD":
                assert protein["sequence"][i] == r
                key = f"{r}{i+1}{a}"
                assert key not in features
                features[key] = list(x[[0,2,9]])
    percentile = {f"{r}{i}{a}":v for g,i,r,a,v,_ in json.loads(SOURCES["endpoint"].read_text()) if g=="G6PD"}
    assert features.keys()==percentile.keys() and len(features)==9676
    with SOURCES["turnover"].open() as fp:
        rows=[r for r in csv.DictReader(fp) if r["endpoint"]=="kcat" and r["added_ligand"]=="not_applicable" and r["status"]=="numeric_available"]
    refs={r["study"]:float(r["value"]) for r in rows if r["background_type"]=="WT"}
    targets=defaultdict(list)
    for r in rows:
        if r["background_type"]=="single":
            targets[r["substitutions"]].append(float(r["value"])/refs[r["study"]])
    keys=sorted(targets)
    positions=lambda k:int(k[1:-1])
    missing = [k for k in keys if k not in features]
    keys=[k for k in keys if k in features]
    weights=Counter(map(positions,keys))
    abundance={}
    for line in SOURCES["reporter"].read_text().splitlines():
        row=json.loads(line)
        if row["position_1based"] in weights:
            k=f"{row['ref_aa']}{row['position_1based']}{row['alt_aa']}"
            if k in percentile and row["abundance_score"] is not None:
                abundance[k]=row["abundance_score"]
    ak=sorted(abundance)
    aw=Counter(map(positions,ak))
    y=[np.median(targets[k]) for k in keys]
    w=[1/weights[positions(k)] for k in keys]
    state=dict(protein_id=protein["protein_id"],reference_sha256=protein["sha256"],
               activity=fit([features[k] for k in keys],y,w),
               activity_control=fit([[percentile[k]] for k in keys],y,w),
               abundance=fit([[percentile[k]] for k in ak],[abundance[k] for k in ak],[1/aw[positions(k)] for k in ak]),
               thresholds=dict(activity=.5,abundance=.75),activity_feature_names=["SaProt_T","ESMC_LLR","ThermoMPNN_ddG"],
               abundance_input="ASPA-fitted relative reporter-abundance percentile fraction",
               training_turnover_variants=keys,training_residues=sorted(weights),training_abundance_rows=len(ak),
               training_missing_features=missing)
    write_new(HERE/"model.json",state)
    output=io.StringIO()
    writer=csv.DictWriter(output,fieldnames=["substitution","training_residue","turnover_fraction","reporter_abundance","phenotype","control_turnover"],lineterminator="\n")
    writer.writeheader()
    for k in sorted(features,key=lambda v:(positions(v),v)):
        writer.writerow(dict(substitution=k,training_residue=int(positions(k) in weights),**predict(state,*features[k],percentile[k])))
    with (HERE/"predictions.csv.gz").open("xb") as fp:
        fp.write(gzip.compress(output.getvalue().encode(),mtime=0))
    write_new(HERE/"score-lock.json", dict(created_utc=datetime.now(timezone.utc).isoformat(), public_predictions=len(features),
              model_sha256=digest(HERE/"model.json"),predictions_sha256=digest(HERE/"predictions.csv.gz"),
              protocol_sha256=digest(HERE/"PROTOCOL.md"),training_code_sha256=digest(Path(__file__)),predictor_sha256=digest(HERE/"predict.py"),
              source_hashes={key:digest(path) for key,path in SOURCES.items()}))

if __name__=="__main__":
    train()
