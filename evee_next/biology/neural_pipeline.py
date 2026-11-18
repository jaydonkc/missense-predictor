"""Attach biological records to an explicitly bound frozen neural readout.

Providers supply already transformed coordinates, their exact feature-map ID,
the resolved genomic/protein identity, and a source role for every head block.
This module extracts no features, loads no artifacts, and fits no models. Source
roles are a caller qualification contract; they cannot prove how an input was
produced. Arithmetic accounting remains development-only and is not mechanism
evidence or a measured mutation effect.
"""

from dataclasses import dataclass, field
import re
from types import MappingProxyType
from typing import Literal, Mapping, Protocol, Sequence

import numpy as np

from evee_next.interpretation.neural import FrozenReluHead, NeuralAccounting, immutable
from evee_next.report import (
    IdentityProvider, MissenseMapping, NormalizedVariant, VariantAnalyzer,
    VariantEvidenceReport, VariantRequest,
)
from .chemistry import Chemistry
from .pipeline import BiologyProvider, ContextRecord, EffectRecord


SourceRole = Literal[
    "DNA", "direct_protein", "amino_acid_chemistry", "availability",
    "reference_structure", "sequence_likelihood", "reference_context_x_chemistry",
]
ClinicalStatus = Literal[
    "invalid_request", "unsupported", "outside_admitted_roster", "mapping_conflict",
    "ambiguous_mapping", "partial", "scored",
]
_INTERPRETATION = {
    "DNA": "task_trained_DNA",
    "direct_protein": "biologically_unresolved_embedding",
    "amino_acid_chemistry": "amino_acid_chemistry",
    "availability": "availability_indicator",
    "reference_structure": "reference_structure_context",
    "sequence_likelihood": "sequence_likelihood_not_mechanism",
    "reference_context_x_chemistry": "reference_context_x_chemistry",
}


@dataclass(frozen=True)
class NeuralFeatureRecord:
    """One identity-bound set of transformed inputs; missing blocks stay absent.

    ``block_sources`` declares every expected block, including unavailable ones.
    Only standalone sources are supported. In particular, an EVEE score, call,
    logit, or feature derived from one must not be declared as another source.
    """

    variant: NormalizedVariant
    mapping: MissenseMapping
    feature_map_id: str
    blocks: Mapping[str, Sequence[float] | np.ndarray]
    block_sources: Mapping[str, SourceRole]
    provenance: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.variant, NormalizedVariant) or not isinstance(self.mapping, MissenseMapping):
            raise ValueError("Neural features require exact normalized genomic and protein identities")
        if not isinstance(self.feature_map_id, str) or re.fullmatch(r"[0-9a-f]{64}", self.feature_map_id) is None:
            raise ValueError("Neural features require a qualified feature-map SHA")
        object.__setattr__(self, "blocks", MappingProxyType({name: immutable(value) for name, value in self.blocks.items()}))
        object.__setattr__(self, "block_sources", MappingProxyType(dict(self.block_sources)))
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


class NeuralFeatureProvider(Protocol):
    def features(self, variant: NormalizedVariant, mapping: MissenseMapping,
                 required_blocks: tuple[str, ...]) -> NeuralFeatureRecord: ...


@dataclass(frozen=True)
class NeuralClinicalReport:
    """The existing unscored identity report plus optional neural accounting.

    ``status`` and ``missing`` describe this scoring adapter. ``identity`` keeps
    the original resolver report, including its unscored identity-only status.
    Provider provenance is separate so it cannot replace reference-source pins.
    """

    identity: VariantEvidenceReport
    status: ClinicalStatus
    score: NeuralAccounting | None
    missing: tuple[str, ...]
    provenance: Mapping[str, str]
    feature_provenance: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class NeuralBiologyReport:
    clinical: NeuralClinicalReport
    chemistry: Mapping
    context: ContextRecord
    effects: tuple[EffectRecord, ...]
    effect_used_in_clinical_score: bool = field(default=False, init=False)
    mechanism_status: Literal["unresolved"] = field(default="unresolved", init=False)


class NeuralBiologyPipeline:
    """Reuse identity resolution and preserve the frozen DNA/protein calculation.

    The feature provider owns extraction, transformation, reference matching and
    qualification. The head and input declarations must contain DNA and direct
    protein pathways. No effect or newly calculated chemistry/context is inserted
    into the score. Existing chemistry/context coordinates in a supplied frozen
    head remain part of that original calculation.
    """

    def __init__(self, identity_provider: IdentityProvider, chemistry: Chemistry, *,
                 head: FrozenReluHead | None = None,
                 feature_provider: NeuralFeatureProvider | None = None,
                 provider: BiologyProvider | None = None):
        if head is not None:
            if not isinstance(head, FrozenReluHead):
                raise ValueError("A frozen ReLU head is required")
            interpretations = {block.interpretation for block in head.blocks}
            if any("evee" in block.name.casefold() for block in head.blocks):
                raise ValueError("An EVEE inference block is not a standalone input")
            if not {"task_trained_DNA", "biologically_unresolved_embedding"} <= interpretations:
                raise ValueError("Clinical model must preserve DNA and direct protein inputs")
        self.analyzer = VariantAnalyzer(identity_provider)
        self.chemistry = chemistry
        self.head = head
        self.feature_provider = feature_provider
        self.provider = provider

    def _clinical(self, identity: VariantEvidenceReport) -> NeuralClinicalReport:
        provenance = dict(identity.provenance)
        if identity.mapping is None or self.head is None:
            return NeuralClinicalReport(identity, identity.status, None, identity.missing, provenance)

        head = self.head
        provenance.update(feature_map_id=head.feature_map_id, head_parameter_sha256=head.parameter_sha256)
        required = tuple(block.name for block in head.blocks)
        if self.feature_provider is None:
            return NeuralClinicalReport(identity, "partial", None, required, provenance)

        record = self.feature_provider.features(identity.normalized_variant, identity.mapping, required)
        if not isinstance(record, NeuralFeatureRecord):
            raise ValueError("A typed identity-bound neural feature record is required")
        if record.variant != identity.normalized_variant or record.mapping != identity.mapping:
            raise ValueError("Neural feature genomic allele or protein mapping differs")
        if record.feature_map_id != head.feature_map_id:
            raise ValueError("Feature coordinate map differs from frozen head")
        if set(record.block_sources) != set(required) or not set(record.blocks) <= set(required):
            raise ValueError("Neural feature block declarations differ from the frozen head")
        for block in head.blocks:
            role = record.block_sources[block.name]
            if role not in _INTERPRETATION or _INTERPRETATION[role] != block.interpretation:
                raise ValueError("Unsupported or inconsistent standalone neural input source")
            if block.name in record.blocks and record.blocks[block.name].shape != (block.end - block.start,):
                raise ValueError(f"Wrong transformed input shape for {block.name}")

        missing = tuple(name for name in required if name not in record.blocks)
        if missing:
            return NeuralClinicalReport(identity, "partial", None, missing, provenance, dict(record.provenance))
        transformed = np.concatenate([record.blocks[name] for name in required])
        score = head.account(transformed, record.feature_map_id)
        return NeuralClinicalReport(identity, "scored", score, (), provenance, dict(record.provenance))

    def analyze(self, request: VariantRequest) -> NeuralBiologyReport:
        identity = self.analyzer.analyze(request)
        clinical = self._clinical(identity)
        if identity.mapping is None:
            return NeuralBiologyReport(clinical, {}, ContextRecord("unavailable"), ())
        mapping = identity.mapping
        chemistry = self.chemistry.calculate(mapping.ref_aa, mapping.alt_aa)
        context = self.provider.context(mapping) if self.provider else ContextRecord("unavailable")
        effects = self.provider.effects(mapping) if self.provider else ()
        if not isinstance(context, ContextRecord) or not isinstance(effects, tuple) or not all(
                isinstance(effect, EffectRecord) for effect in effects):
            raise ValueError("Biological providers must return typed context and effect records")
        return NeuralBiologyReport(clinical, chemistry, context, effects)
