"""Typed missense evidence reports and direct-readout score accounting."""

from dataclasses import dataclass, field
import math
import re
from typing import Literal, Mapping, Protocol, Sequence

import numpy as np
from scipy.special import expit


AA = "ACDEFGHIKLMNPQRSTVWY"
CHROMOSOMES = {str(n) for n in range(1, 23)} | {"X", "Y"}


@dataclass(frozen=True)
class VariantRequest:
    chromosome: str
    position_1based: int
    ref: str
    alt: str
    transcript_id: str | None = None
    protein_id: str | None = None
    organism: str = "human"
    assembly: str = "GRCh38"


@dataclass(frozen=True)
class NormalizedVariant:
    chromosome: str
    position_1based: int
    position_0based: int
    ref: str
    alt: str
    variant_id: str
    organism: str = "Homo sapiens"
    assembly: str = "GRCh38"


@dataclass(frozen=True)
class MissenseMapping:
    transcript_id: str
    protein_id: str
    strand: Literal["+", "-"]
    residue_1based: int
    ref_aa: str
    alt_aa: str

    def __post_init__(self):
        if (not re.fullmatch(r"ENST\d+\.\d+", self.transcript_id)
                or not re.fullmatch(r"ENSP\d+\.\d+", self.protein_id)):
            raise ValueError("The missense route requires versioned Ensembl mappings")
        if (self.strand not in {"+", "-"} or type(self.residue_1based) is not int or self.residue_1based < 1
                or len(self.ref_aa) != 1 or self.ref_aa not in AA
                or len(self.alt_aa) != 1 or self.alt_aa not in AA or self.ref_aa == self.alt_aa):
            raise ValueError("The mapping must identify a canonical single-amino-acid substitution")


@dataclass(frozen=True)
class EvidenceRecord:
    kind: Literal["context", "measured_effect", "predicted_effect", "reported_effect", "hypothesis"]
    statement: str
    source: str


class IdentityProvider(Protocol):
    provenance: Mapping[str, str]

    def mappings(self, variant: NormalizedVariant) -> tuple[MissenseMapping, ...]: ...


@dataclass(frozen=True)
class IdentityResolution:
    """Optional richer response from an identity provider."""

    mappings: tuple[MissenseMapping, ...]
    applicability: str
    empty_status: Literal["invalid_request", "unsupported", "outside_admitted_roster", "mapping_conflict", "partial"] = "partial"
    missing: tuple[str, ...] = ()
    conflicts: tuple[str, ...] = ()
    unresolved_questions: tuple[str, ...] = ()
    provenance: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class FeatureRecord:
    blocks: Mapping[str, Sequence[float] | np.ndarray]
    provenance: Mapping[str, str] = field(default_factory=dict)


class FeatureProvider(Protocol):
    def features(self, variant: NormalizedVariant, mapping: MissenseMapping,
                 required_blocks: tuple[str, ...]) -> FeatureRecord: ...


@dataclass(frozen=True)
class BlockContribution:
    block: str
    signed_logit: float
    interpretation: Literal["biologically_unresolved_embedding", "substitution_chemistry", "archived_evee_logit"]


@dataclass(frozen=True)
class PathogenicityScore:
    value: float
    logit: float
    intercept: float
    contributions: tuple[BlockContribution, ...]
    cutoff: float
    predicted_class: Literal["pathogenic", "benign"]
    model_arm: str
    calibration_status: Literal["uncalibrated_sigmoid_score"] = "uncalibrated_sigmoid_score"
    evaluation_scope: Literal["development_only"] = "development_only"


@dataclass(frozen=True)
class VariantEvidenceReport:
    request: VariantRequest
    status: Literal["invalid_request", "unsupported", "outside_admitted_roster", "mapping_conflict",
                    "ambiguous_mapping", "partial", "scored"]
    normalized_variant: NormalizedVariant | None = None
    mapping: MissenseMapping | None = None
    mapping_options: tuple[MissenseMapping, ...] = ()
    score: PathogenicityScore | None = None
    evidence: tuple[EvidenceRecord, ...] = ()
    missing: tuple[str, ...] = ()
    conflicts: tuple[str, ...] = ()
    unresolved_questions: tuple[str, ...] = ()
    provenance: Mapping[str, str] = field(default_factory=dict)
    applicability: str = "outside_scoring_scope"
    mechanism_status: Literal["unresolved"] = "unresolved"


@dataclass(frozen=True)
class FrozenDirectReadout:
    arm: str
    blocks: tuple[tuple[str, int], ...]
    means: np.ndarray
    scales: np.ndarray
    weights: np.ndarray
    intercept: float
    cutoff: float
    provenance: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self):
        names = {name for name, _ in self.blocks}
        if (not self.blocks or len(names) != len(self.blocks) or not names <= {"D", "P", "ESMC", "C", "EVEE"}
                or any(type(width) is not int or width < 1 for _, width in self.blocks)):
            raise ValueError("Invalid direct-model feature blocks")
        size = sum(width for _, width in self.blocks)
        for name in ("means", "scales", "weights"):
            array = np.array(getattr(self, name), dtype=np.float64, copy=True)
            if array.shape != (size,) or not np.isfinite(array).all():
                raise ValueError(f"Invalid model {name}")
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if "EVEE" in names and (self.blocks[0] != ("EVEE", 1) or self.means[0] != 0
                               or self.scales[0] != 1 or self.weights[0] != 1 or self.intercept != 0):
            raise ValueError("Archived EVEE must be a fixed logit offset")
        if (np.any(self.scales <= 0) or not math.isfinite(self.intercept)
                or not 0 <= self.cutoff <= math.nextafter(1.0, math.inf)):
            raise ValueError("Invalid model scale, intercept, or cutoff")

    def score(self, features: Mapping[str, Sequence[float] | np.ndarray]) -> PathogenicityScore:
        offset = 0
        contributions = []
        for block, width in self.blocks:
            vector = np.asarray(features[block], dtype=np.float64)
            if vector.shape != (width,) or not np.isfinite(vector).all():
                raise ValueError(f"Invalid feature block {block}")
            section = slice(offset, offset + width)
            offset += width
            terms = ((vector - self.means[section]) / self.scales[section]) * self.weights[section]
            logit = float(terms.sum())
            if not math.isfinite(logit):
                raise ValueError(f"Nonfinite contribution from feature block {block}")
            kind = {"C": "substitution_chemistry", "EVEE": "archived_evee_logit"}.get(
                block, "biologically_unresolved_embedding")
            contributions.append(BlockContribution(block, logit, kind))
        logit = self.intercept + sum(term.signed_logit for term in contributions)
        if not math.isfinite(logit):
            raise ValueError("Nonfinite combined logit")
        value = float(expit(logit))
        return PathogenicityScore(value, logit, self.intercept, tuple(contributions), self.cutoff,
                                  "pathogenic" if value >= self.cutoff else "benign", self.arm)


def normalize(request: VariantRequest) -> NormalizedVariant:
    if request.organism.strip().casefold() not in {"human", "homo sapiens", "9606"}:
        raise NotImplementedError("This route supports human variants")
    if request.assembly.strip().casefold() != "grch38":
        raise NotImplementedError("This route requires the GRCh38 assembly")
    chromosome = request.chromosome.strip().upper().removeprefix("CHR")
    if chromosome not in CHROMOSOMES:
        raise NotImplementedError("This route supports canonical nuclear GRCh38 chromosomes")
    if type(request.position_1based) is not int or request.position_1based < 1:
        raise ValueError("Use a positive, one-based genomic position")
    ref, alt = request.ref.strip().upper(), request.alt.strip().upper()
    if len(ref) != 1 or len(alt) != 1:
        raise NotImplementedError("This route supports single-nucleotide substitutions")
    if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        raise ValueError("Use distinct canonical reference and alternate bases")
    index = request.position_1based - 1
    return NormalizedVariant(chromosome, request.position_1based, index, ref, alt,
                             f"chr{chromosome}:{index}:{ref}:{alt}")


def matching_mappings(request: VariantRequest, options: tuple[MissenseMapping, ...]) -> tuple[MissenseMapping, ...]:
    transcript = request.transcript_id.strip().upper() if request.transcript_id is not None else None
    protein = request.protein_id.strip().upper() if request.protein_id is not None else None
    return tuple(row for row in options if transcript in {None, row.transcript_id} and protein in {None, row.protein_id})


@dataclass(frozen=True)
class VariantAnalyzer:
    identity_provider: IdentityProvider
    feature_provider: FeatureProvider | None = None
    model: FrozenDirectReadout | None = None
    model_provenance: Mapping[str, str] = field(default_factory=dict)

    def analyze(self, request: VariantRequest) -> VariantEvidenceReport:
        provenance = dict(self.identity_provider.provenance)
        try:
            variant = normalize(request)
        except NotImplementedError as exc:
            return VariantEvidenceReport(request, "unsupported", unresolved_questions=(str(exc),), provenance=provenance)
        except ValueError as exc:
            return VariantEvidenceReport(request, "invalid_request", conflicts=(str(exc),), provenance=provenance)

        detailed = getattr(self.identity_provider, "resolve_identity", None)
        resolution = detailed(variant) if callable(detailed) else None
        options = resolution.mappings if resolution is not None else self.identity_provider.mappings(variant)
        if resolution is not None:
            provenance.update(resolution.provenance)
        if not options:
            if resolution is not None:
                return VariantEvidenceReport(request, resolution.empty_status, variant, missing=resolution.missing,
                                             conflicts=resolution.conflicts, unresolved_questions=resolution.unresolved_questions,
                                             provenance=provenance, applicability=resolution.applicability)
            return VariantEvidenceReport(request, "outside_admitted_roster", variant, missing=("admitted_identity",),
                                         unresolved_questions=("A verified missense mapping and feature cache are required.",),
                                         provenance=provenance)
        matches = matching_mappings(request, options)
        if not matches:
            return VariantEvidenceReport(request, "mapping_conflict", variant, mapping_options=options,
                                         conflicts=("The requested transcript/protein versions do not match the verified options.",),
                                         provenance=provenance)
        if len(matches) > 1:
            return VariantEvidenceReport(request, "ambiguous_mapping", variant, mapping_options=matches,
                                         missing=("resolved_transcript_and_protein",), provenance=provenance)

        mapping = matches[0]
        provenance.update(self.model_provenance)
        evidence = EvidenceRecord("context", f"{mapping.protein_id} carries the mapped substitution "
                                  f"{mapping.ref_aa}{mapping.residue_1based}{mapping.alt_aa}.",
                                  provenance.get("identity_mapping", "injected_identity_provider"))
        report = dict(normalized_variant=variant, mapping=mapping, evidence=(evidence,), provenance=provenance,
                      applicability=resolution.applicability if resolution is not None else "admitted_development_missense",
                      unresolved_questions=("What biological effect does this substitution cause in the relevant context?",))
        if self.model is None:
            return VariantEvidenceReport(request, "partial", missing=("selected_development_model",), **report)
        provenance.update(self.model.provenance)
        required = tuple(name for name, _ in self.model.blocks)
        if self.feature_provider is None:
            return VariantEvidenceReport(request, "partial", missing=required, **report)
        feature_record = self.feature_provider.features(variant, mapping, required)
        provenance.update(feature_record.provenance)
        missing = tuple(name for name in required if name not in feature_record.blocks)
        if missing:
            return VariantEvidenceReport(request, "partial", missing=missing, **report)
        return VariantEvidenceReport(request, "scored", score=self.model.score(feature_record.blocks), **report)
