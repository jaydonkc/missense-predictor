"""Outcome-free identities, evidence contracts and public report serialization."""
from dataclasses import dataclass, field, fields, is_dataclass
import hashlib
import json
import math
import re
from types import MappingProxyType
from typing import Mapping, Protocol

from evee_next.report import MissenseMapping, NormalizedVariant, VariantRequest, matching_mappings


def canonical_bytes(value):
    return json.dumps(public_value(value), sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def canonical_sha(value):
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def require_sha(value):
    if not isinstance(value, str) or re.fullmatch('[a-f0-9]{64}', value) is None:
        raise ValueError('An exact SHA256 is required')


def _frozen(value):
    if isinstance(value, Mapping):
        return MappingProxyType({str(k): _frozen(v) for k, v in value.items()})
    if isinstance(value, (tuple, list)): return tuple(_frozen(v) for v in value)
    if value is None or isinstance(value, (str, bool, int)): return value
    if isinstance(value, float) and math.isfinite(value): return value
    raise ValueError('Only finite JSON values belong in evidence metadata')


def public_value(value):
    if isinstance(value, BaselineState) or isinstance(value, BaselineResult):
        raise ValueError('Private baseline state cannot be serialized')
    if is_dataclass(value): return {f.name: public_value(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, Mapping): return {k: public_value(v) for k,v in value.items()}
    if isinstance(value, (tuple, list)): return [public_value(v) for v in value]
    return value


@dataclass(frozen=True)
class VariantContext:
    status: str
    detail: str
    variant: NormalizedVariant | None
    mappings: tuple[MissenseMapping, ...] = ()
    mapping: MissenseMapping | None = None
    transcripts: tuple = ()
    provenance: Mapping[str, str] = field(default_factory=dict)
    genomic_validated: bool = False

    def __post_init__(self):
        object.__setattr__(self, 'mappings', tuple(self.mappings))
        object.__setattr__(self, 'transcripts', tuple(self.transcripts))
        object.__setattr__(self, 'provenance', _frozen(self.provenance))
        if self.mapping is not None and (self.mapping not in self.mappings or len(self.mappings) != 1):
            raise ValueError('A selected mapping must be unique and identity-bound')
        if self.genomic_validated and self.variant is None:
            raise ValueError('Validated context requires a normalized allele')

    @classmethod
    def from_request(cls, request, resolver):
        resolution = resolver.resolve(request)
        variant = resolution.variant
        options = matching_mappings(request, resolution.mappings)
        status, detail = resolution.status, resolution.detail
        if resolution.mappings and not options:
            status, detail = 'mapping_conflict', 'Requested transcript/protein does not match the reference'
        elif len(options) > 1:
            status = 'ambiguous_mapping'
        elif not options and any(a.status in {'synonymous','stop_gained','stop_lost','stop_retained',
                        'initiation_context_requires_review'} for a in resolution.assessments):
            status = 'unsupported_consequence'
        mapping = options[0] if len(options) == 1 else None
        bad = {'invalid_request', 'reference_mismatch', 'reference_integrity_error', 'missing_reference',
               'unsupported', 'unsupported_consequence', 'mapping_conflict'}
        valid = variant is not None and status not in bad
        transcripts = resolver.transcript_context(variant) if valid and hasattr(resolver, 'transcript_context') else ()
        return cls(status, detail, variant, options, mapping, transcripts, resolution.provenance, valid)


@dataclass(frozen=True)
class ExecutionPolicy:
    mode: str = 'cache-only'
    per_variant_bytes: int = 256 * 1024**2
    invocation_bytes: int = 1024**3
    task_seconds: float = 300.
    invocation_seconds: float = 900.
    transport_retries: int = 2
    allowed_runner_ids: tuple[str, ...] = ()
    allowed_hosts: tuple[str, ...] = ()
    authorization_refs: tuple[str, ...] = ()

    def __post_init__(self):
        if self.mode not in {'cache-only', 'bounded-retrieval'}:
            raise ValueError('Unknown execution policy')
        for name in ('per_variant_bytes', 'invocation_bytes'):
            if type(getattr(self,name)) is not int or getattr(self,name) <= 0:
                raise ValueError('Positive integer byte budgets required')
        for name in ('task_seconds', 'invocation_seconds'):
            x = getattr(self, name)
            if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x <= 0:
                raise ValueError('Finite positive time budgets required')
        if type(self.transport_retries) is not int or not 0 <= self.transport_retries <= 2:
            raise ValueError('At most two transport retries are supported')
        for name in ('allowed_runner_ids', 'allowed_hosts', 'authorization_refs'):
            vals = tuple(getattr(self,name))
            if any(not isinstance(v,str) or not v for v in vals) or len(set(vals)) != len(vals):
                raise ValueError('Unique nonempty execution declarations required')
            object.__setattr__(self,name,vals)


_IDENTITY_KEYS = frozenset({'assembly','organism','chromosome','position_1based','ref','alt',
    'transcript_id','protein_id','protein_sha256','residue_1based','ref_aa','alt_aa','strand',
    'gene_id','gene_symbol','reference_sha256','annotation_sha256','fragment_id','chain_id'})
_FORBIDDEN = frozenset({'clinical_label','label','target','targets','known_error','evee_score','evee_call'})


def _check_metadata(value):
    if isinstance(value, Mapping):
        if _FORBIDDEN.intersection(str(k).lower() for k in value):
            raise ValueError('Outcome or comparator inputs are not allowed in inference requests')
        for v in value.values(): _check_metadata(v)
    elif isinstance(value,(tuple,list)):
        for v in value: _check_metadata(v)


@dataclass(frozen=True)
class EvidenceRequest:
    provider_id: str
    recipe_sha256: str
    identity: Mapping
    parameters: Mapping = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.provider_id,str) or re.fullmatch('[a-z][a-z0-9_.-]*',self.provider_id) is None:
            raise ValueError('Stable provider ID required')
        require_sha(self.recipe_sha256)
        if not self.identity or not set(self.identity) <= _IDENTITY_KEYS:
            raise ValueError('Unsupported or empty evidence identity')
        _check_metadata(self.parameters)
        object.__setattr__(self,'identity',_frozen(self.identity))
        object.__setattr__(self,'parameters',_frozen(self.parameters))

    @property
    def key(self): return canonical_sha(self)


_STATUSES = frozenset({'available','partial','not_applicable','unavailable','ambiguous','failed','pending_authorization'})


@dataclass(frozen=True)
class EvidenceDecision:
    request: EvidenceRequest
    applicable: bool
    status: str
    reason: str
    values: Mapping[str, float | int | None] = field(default_factory=dict)
    units: Mapping[str, str] = field(default_factory=dict)
    provenance: Mapping[str, str] = field(default_factory=dict)
    metadata: Mapping = field(default_factory=dict)
    uncertainty: str = 'not_quantified'
    qualification: str = 'diagnostic'
    cache_disposition: str = 'not_used'
    execution: str = 'not_run'
    scoring_used: tuple[str, ...] = ()

    def __post_init__(self):
        if not isinstance(self.request,EvidenceRequest) or type(self.applicable) is not bool:
            raise ValueError('Typed evidence request and applicability required')
        if self.status not in _STATUSES or not isinstance(self.reason,str) or not self.reason:
            raise ValueError('Explicit evidence status and reason required')
        if set(self.units) != set(self.values) or any(not isinstance(u,str) or not u for u in self.units.values()):
            raise ValueError('Every evidence field needs an explicit unit')
        for x in self.values.values():
            if x is not None and (isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x)):
                raise ValueError('Evidence must be finite numerical values or explicit null')
        known = sum(x is not None for x in self.values.values())
        if self.status == 'available' and (not known or known != len(self.values) or not self.provenance):
            raise ValueError('Available evidence requires complete finite values and provenance')
        if self.status == 'partial' and (not known or not self.provenance):
            raise ValueError('Partial evidence needs some qualified values and provenance')
        if self.status not in {'available','partial'} and known:
            raise ValueError('Unusable evidence cannot contain usable numerical values')
        if self.status == 'not_applicable' and self.applicable:
            raise ValueError('Not-applicable evidence cannot be applicable')
        if not set(self.scoring_used) <= set(self.values) or any(self.values[k] is None for k in self.scoring_used):
            raise ValueError('Only available evidence fields can be consumed')
        for name in ('values','units','provenance','metadata'):
            object.__setattr__(self,name,_frozen(getattr(self,name)))
        object.__setattr__(self,'scoring_used',tuple(self.scoring_used))


@dataclass(frozen=True)
class BaselinePrediction:
    status: str
    score: float | None = None
    cutoff: float | None = None
    checkpoint: str | None = None
    feature_map_id: str | None = None
    scope: str = 'unscored'
    score_type: str = 'uncalibrated_sigmoid'
    reason: str = ''

    def __post_init__(self):
        if self.score is not None:
            if isinstance(self.score,bool) or not math.isfinite(self.score) or not 0 <= self.score <= 1:
                raise ValueError('Finite 0..1 score required')
            if self.cutoff is None or isinstance(self.cutoff,bool) or not math.isfinite(self.cutoff) or not 0 <= self.cutoff <= 1:
                raise ValueError('A fixed operating point is required for a score')
            if not self.checkpoint or not self.feature_map_id: raise ValueError('Score needs checkpoint/feature identity')
        elif self.status == 'scored': raise ValueError('Scored status needs a score')

    @property
    def call(self): return None if self.score is None else self.score >= self.cutoff


@dataclass(frozen=True)
class BaselineState:
    variant_key: str
    checkpoint: str
    feature_map_id: str
    a0: tuple[float, ...]
    W1: tuple[float, ...]
    b1: float

    def __post_init__(self):
        if not self.a0 or len(self.a0) != len(self.W1) or not all(math.isfinite(x) for x in (*self.a0,*self.W1,self.b1)):
            raise ValueError('Finite identity-bound baseline state required')
        object.__setattr__(self,'a0',tuple(self.a0)); object.__setattr__(self,'W1',tuple(self.W1))


@dataclass(frozen=True)
class BaselineResult:
    prediction: BaselinePrediction
    state: BaselineState | None = None

    def __post_init__(self):
        if self.state is not None and (self.prediction.checkpoint,self.prediction.feature_map_id) != (self.state.checkpoint,self.state.feature_map_id):
            raise ValueError('Private state differs from public checkpoint')


@dataclass(frozen=True)
class RoutedBiologyReport:
    context: VariantContext
    baseline: BaselinePrediction
    final: BaselinePrediction
    evidence: tuple[EvidenceDecision, ...]
    fusion_status: str = 'not_configured'
    mechanism_status: str = 'unresolved'
    schema: str = 'routed-missense-biology/v1'


def report_dict(report):
    result = public_value(report)
    result['baseline']['call'] = report.baseline.call
    result['final']['call'] = report.final.call
    return result


def validate_request_identity(request,context,*,selected=False):
    """Check an evidence request against the allele being reported/scored."""
    v=context.variant
    if v is None:raise ValueError('Evidence lacks a normalized genomic context')
    identity=request.identity
    for key in ('assembly','chromosome','position_1based','ref','alt'):
        if identity.get(key)!=getattr(v,key):raise ValueError('Evidence differs from current genomic '+key)
    if identity.get('organism',v.organism)!=v.organism:raise ValueError('Evidence organism differs')
    mappings=(context.mapping,) if selected else context.mappings
    mappings=tuple(m for m in mappings if m is not None)
    protein_keys=('transcript_id','protein_id','strand','residue_1based','ref_aa','alt_aa')
    if any(k in identity for k in protein_keys):
        matches=[m for m in mappings if all(identity.get(k,getattr(m,k))==getattr(m,k) for k in protein_keys)]
        # Transcript-only RNA evidence can remain valid with no protein mapping.
        if not matches:
            tx=identity.get('transcript_id')
            transcript_matches=[t for t in context.transcripts if t.transcript_id==tx
                and ('strand' not in identity or identity['strand']==t.strand)]
            if selected or not transcript_matches or any(k in identity for k in ('residue_1based','ref_aa','alt_aa')):
                raise ValueError('Evidence differs from selected protein/transcript context')
    if 'protein_sha256' in identity:
        refs=[t for t in context.transcripts if getattr(t,'protein_id',None)==identity.get('protein_id')]
        hashes={getattr(t,'reference_hashes',{}).get('protein_sha256') for t in refs}
        hashes.update(hashlib.sha256(t.protein_sequence.encode()).hexdigest() for t in refs if getattr(t,'protein_sequence',None))
        if identity['protein_sha256'] not in hashes:raise ValueError('Evidence protein reference differs')
    if 'reference_sha256' in identity:
        known=set(context.provenance.values())
        for t in context.transcripts:known.update(getattr(t,'reference_hashes',{}).values())
        if identity['reference_sha256'] not in known:raise ValueError('Evidence genomic reference is unqualified')


class EvidenceProvider(Protocol):
    provider_id: str
    def plan(self, context: VariantContext) -> tuple[EvidenceRequest, ...]: ...
    def obtain(self, request: EvidenceRequest, store, policy: ExecutionPolicy) -> EvidenceDecision: ...


class BaselineProvider(Protocol):
    def predict(self, context: VariantContext) -> BaselineResult: ...
