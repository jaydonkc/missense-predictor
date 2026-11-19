"""Score one prepared GRCh38 missense substitution with frozen v1 heads."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
from scipy.special import expit, logit

ASSETS = Path(__file__).resolve().parent / "assets"
BASES = "ACGT"
AA = "ACDEFGHIKLMNPQRSTVWY"
PROFILE = (
    "effective_species", "weighted_ref_fraction", "weighted_alt_fraction", "entropy_nats",
    "log_ref_alt", "blosum_log_ref_alt", "nearest_callable_distance", "nearest_alt_support",
    "clades_callable", "clades_with_alt", "primate_callable", "primate_alt",
)
COORDINATES = {"assembly", "chrom", "position_0based", "ref", "alt"}
COMP = str.maketrans("ACGT", "TGCA")


def checked(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"Model asset failed its SHA-256 check: {path.name}")
    return raw


def number(value):
    if value is None:
        return float("nan")
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError("Expected a finite number or explicit null")
    return float(value)


def one_hot(value, alphabet):
    if value not in alphabet:
        return [np.nan] * len(alphabet)
    return [float(char == value) for char in alphabet]


def dna_features(variant):
    ref, alt = variant["ref"], variant["alt"]
    window = variant["dna_window"]
    if window is None:
        window = "N" * 11
    elif not isinstance(window, str) or not re.fullmatch("[ACGTN]{11}", window) or window[5] != ref:
        raise ValueError("dna_window must be null or an 11-base GRCh38 window centered on ref")
    if ref in "AG":
        ref, alt = ref.translate(COMP), alt.translate(COMP)
        window = window.translate(COMP)[::-1]
    flank = window[:5] + window[6:]
    gc = sum(x in "GC" for x in flank) / 10 if "N" not in flank else np.nan
    cpg = int("CG" in window[4:7]) if "N" not in window[4:7] else np.nan
    return [number(variant["phylop"]), *one_hot(ref, BASES), *one_hot(alt, BASES), gc, cpg,
            *(x for base in flank for x in one_hot(base, BASES))]


def protein_features(protein, chemistry):
    required = {"ref_aa", "alt_aa", "transcript_id", "protein_id", "residue_1based"}
    if not isinstance(protein, dict) or set(protein) not in (required, required | {"alignment"}):
        raise ValueError("Invalid protein feature schema")
    ref, alt = protein["ref_aa"], protein["alt_aa"]
    if ref not in set(AA) or alt not in set(AA) or ref == alt:
        raise ValueError("Requested protein substitution must change a canonical amino acid")
    if type(protein["residue_1based"]) is not int or protein["residue_1based"] <= 0:
        raise ValueError("Invalid protein residue number")
    for key, pattern in (("protein_id", r"ENSP[0-9]+\.[0-9]+"), ("transcript_id", r"ENST[0-9]+\.[0-9]+")):
        if not isinstance(protein[key], str) or not re.fullmatch(pattern, protein[key]):
            raise ValueError(f"{key} must be a versioned Ensembl identifier")
    deltas = [chemistry[k][alt] - chemistry[k][ref] for k in ("charge", "volume", "polarity", "hydrophobicity")]
    features = [chemistry["grantham"][ref][alt], chemistry["pam250"][ref][alt], *deltas, *map(abs, deltas),
                *one_hot(ref, AA), *one_hot(alt, AA)]
    profile = protein.get("alignment")
    if profile is None:
        return features + [np.nan] * len(PROFILE), "DNA and protein chemistry"
    if not isinstance(profile, dict) or set(profile) != set(PROFILE):
        raise ValueError("Expected the twelve original allele-specific profile statistics")
    values = [number(profile[k]) for k in PROFILE]
    if not np.all(np.isfinite(values)):
        raise ValueError("A qualified alignment profile must have twelve measured feature values")
    return features + values, "DNA and aligned protein"


def parse_heads(raw):
    heads = []
    for head in raw:
        trees = [tuple(np.asarray(column) for column in tree) for tree in head["trees"]]
        heads.append((head["base"], trees))
    return heads


def raw_margin(head, x):
    base, trees = head
    scores = np.full(len(x), base, dtype=float)
    for feature, threshold, miss_left, left, right, value, is_leaf in trees:
        cursor = np.zeros(len(x), dtype=np.int32)
        while True:
            rows = np.flatnonzero(~is_leaf[cursor].astype(bool))
            if not len(rows):
                break
            nodes = cursor[rows]
            v = x[rows, feature[nodes]]
            go_left = np.where(np.isnan(v), miss_left[nodes], v <= threshold[nodes])
            cursor[rows] = np.where(go_left, left[nodes], right[nodes])
        scores += value[cursor]
    return scores


class Predictor:
    def __init__(self):
        self.receipt = json.loads((ASSETS / "manifest.json").read_bytes())
        state = json.loads(gzip.decompress(checked(ASSETS / "model.json.gz", self.receipt["model_sha256"])))
        self.chemistry = json.loads(checked(ASSETS / "chemistry.json", self.receipt["chemistry_sha256"]))
        self.dna = parse_heads(state["DNA only"])
        self.protein = parse_heads(state["protein"])
        self.context = state["full_context"]
        self.full = parse_heads(self.context.pop("heads"))

    def full_context_scores(self, fold, parents, features):
        """Score qualified rows for one frozen head, in manifest feature order."""
        parents = np.asarray(parents, dtype=float)
        features = np.asarray(features, dtype=float)
        if type(fold) is not int or fold not in range(3):
            raise ValueError("Expected a frozen fold from 0 to 2")
        if parents.ndim != 1 or features.shape != (parents.size, len(self.context["features"])):
            raise ValueError("Full-context feature matrix must align with parent scores")
        if not np.isfinite(parents).all() or np.any((parents < 0) | (parents > 1)):
            raise ValueError("Expected finite parent scores between zero and one")
        if np.isinf(features).any():
            raise ValueError("Use NaN for qualified missing features")
        base = logit(np.clip(parents, 1e-12, 1 - 1e-12)) - logit(self.context["parent_cutoff"])
        margin = base + raw_margin(self.full[fold], features) - .5 * self.context["calibration"][fold]
        return expit(margin / 2)

    def score(self, variant):
        if not isinstance(variant, dict) or not COORDINATES <= set(variant) or variant["assembly"] != "GRCh38":
            raise ValueError("Expected a prepared GRCh38 variant")
        chrom, pos, ref, alt = (variant[k] for k in ("chrom", "position_0based", "ref", "alt"))
        if not isinstance(chrom, str) or not re.fullmatch(r"chr([1-9]|1[0-9]|2[0-2]|[XYM])", chrom):
            raise ValueError("Expected a primary GRCh38 chromosome")
        if type(pos) is not int or pos < 0 or ref not in set(BASES) or alt not in set(BASES) or ref == alt:
            raise ValueError("Expected one literal substitution at a nonnegative zero-based coordinate")
        if "full_context" in variant:
            if set(variant) != COORDINATES | {"full_context"}:
                raise ValueError("Use just the qualified full_context bundle with this route")
            context = variant["full_context"]
            if not isinstance(context, dict) or set(context) != {"held_fold", "parent_score", "features"}:
                raise ValueError("Invalid full-context bundle")
            fold, parent = context["held_fold"], number(context["parent_score"])
            if type(fold) is not int or fold not in range(3) or not 0 <= parent <= 1:
                raise ValueError("Expected a qualified held fold and finite frozen parent score")
            cols = self.context["features"]
            if set(context["features"]) != set(cols):
                raise ValueError("Full-context fields must exactly match the frozen feature schema")
            x = np.asarray([[number(context["features"][key]) for key in cols]], dtype=float)
            score = float(self.full_context_scores(fold, np.array([parent]), x)[0])
            route = "original full context"
        else:
            required = COORDINATES | {"phylop", "dna_window"}
            if set(variant) not in (required, required | {"protein"}):
                raise ValueError("Expected dna_window and phylop plus optional qualified protein features")
            values = dna_features(variant)
            heads, route = self.dna, "DNA only"
            if variant.get("protein") is not None:
                additional, route = protein_features(variant["protein"], self.chemistry)
                values.extend(additional)
                heads = self.protein
            x = np.asarray([values], dtype=np.float32)
            score = float(np.mean([expit(raw_margin(head, x)[0]) for head in heads]))
        cutoff = self.receipt["cutoff"]
        return {"variant_id": f"{chrom}:{pos}:{ref}:{alt}", "route": route, "score": score,
                "cutoff": cutoff, "screen_positive": score >= cutoff,
                "model_sha256": self.receipt["model_sha256"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    args = parser.parse_args()
    try:
        report = Predictor().score(json.loads(args.request.read_bytes()))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, str(exc) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
