"""Compose existing genomic resolution/scoring with explicit biological stages."""
from dataclasses import dataclass, field
import math
from typing import Mapping, Protocol
from evee_next.report import VariantAnalyzer, VariantEvidenceReport, VariantRequest, MissenseMapping
from .chemistry import Chemistry


@dataclass(frozen=True)
class ContextRecord:
    status: str
    values: Mapping[str, float | None] = field(default_factory=dict)
    sources: Mapping[str, str] = field(default_factory=dict)
    interpretation: str = 'reference_context_not_mutation_effect'

    def __post_init__(self):
        if self.status not in {'available','partial','unavailable'}:
            raise ValueError('Unknown context status')
        if self.status=='unavailable' and self.values:
            raise ValueError('Unavailable context cannot contain numerical values')
        if self.status!='unavailable' and (not self.values or not self.sources):
            raise ValueError('Available context requires values and provenance')
        known=sum(x is not None for x in self.values.values())
        if self.status=='available' and known!=len(self.values):
            raise ValueError('Available context cannot contain unknown values')
        if self.status=='partial' and not 0<known<len(self.values):
            raise ValueError('Partial context requires known and unknown values')
        if any(x is not None and (isinstance(x,bool) or not math.isfinite(x)) for x in self.values.values()):
            raise ValueError('Nonfinite context')


@dataclass(frozen=True)
class EffectRecord:
    endpoint: str
    value: float
    units: str
    source: str
    qualification: str
    uncertainty: str
    clinical_meaning: str = 'not_a_pathogenicity_probability'

    def __post_init__(self):
        if (not all((self.endpoint,self.units,self.source,self.qualification,self.uncertainty))
                or isinstance(self.value,bool) or not math.isfinite(self.value)):
            raise ValueError('Finite named effect and explicit provenance/qualification required')


class BiologyProvider(Protocol):
    def context(self,mapping: MissenseMapping) -> ContextRecord: ...
    def effects(self,mapping: MissenseMapping) -> tuple[EffectRecord,...]: ...


@dataclass(frozen=True)
class BiologyReport:
    clinical: VariantEvidenceReport
    chemistry: Mapping
    context: ContextRecord
    effects: tuple[EffectRecord,...]
    effect_used_in_clinical_score: bool = False
    mechanism_status: str = 'unresolved'


class BiologyPipeline:
    """Keep score and molecular evidence distinct until integration is qualified.

    An existing standalone analyzer owns nucleotide/transcript/codon validation.
    No effect output silently modifies its clinical score.
    """
    def __init__(self,analyzer: VariantAnalyzer,chemistry: Chemistry,provider: BiologyProvider | None=None):
        model=analyzer.model
        if model is not None:
            names={x[0] for x in model.blocks}
            if 'EVEE' in names:
                raise ValueError('EVEE-dependent clinical model is not standalone')
            if not {'D','ESMC'}<=names and not {'D','P'}<=names:
                raise ValueError('Clinical model must preserve DNA and direct protein inputs')
        self.analyzer=analyzer;self.chemistry=chemistry;self.provider=provider

    def analyze(self,request: VariantRequest) -> BiologyReport:
        clinical=self.analyzer.analyze(request)
        if clinical.mapping is None:
            return BiologyReport(clinical,{},ContextRecord('unavailable'),())
        m=clinical.mapping
        chemistry=self.chemistry.calculate(m.ref_aa,m.alt_aa)
        context=self.provider.context(m) if self.provider else ContextRecord('unavailable')
        effects=self.provider.effects(m) if self.provider else ()
        return BiologyReport(clinical,chemistry,context,effects)
