"""Offline, admitted-input evaluation; this module never discovers targets.

Native assay measurements retain their own units and declared target meaning.
These functions measure already frozen predictions; they cannot fit, select an
assessment threshold, release a cohort, or confer registry qualification.
"""
from dataclasses import dataclass, field
from itertools import combinations
import math
from typing import Mapping

import numpy as np

from .routing import _frozen, require_sha


def _number(value):
    return not isinstance(value,bool) and isinstance(value,(int,float)) and math.isfinite(value)


def _scores(value, *, probability=False):
    if any(not isinstance(k,str) or not k or (v is not None and
           (not _number(v) or (probability and not 0 <= v <= 1))) for k,v in value.items()):
        raise ValueError('Named finite predictions or explicit missing values required')
    return _frozen(value)


@dataclass(frozen=True)
class ComponentContract:
    """Frozen row/head/component and arm-lineage admission, supplied externally."""
    role_map_sha256: str
    assessment_components: Mapping[str,str]
    assessment_heads: Mapping[str,str]
    fitting_components_by_head: Mapping[str,tuple[str,...]]
    prediction_lineage_sha256: Mapping[str,str]
    frozen_receipt_sha256: str
    assessment_source_sha256: str
    arm_interpretations: Mapping[str,str] = field(default_factory=dict)
    prediction_directions: Mapping[str,int] = field(default_factory=dict)

    def __post_init__(self):
        require_sha(self.role_map_sha256); require_sha(self.frozen_receipt_sha256)
        require_sha(self.assessment_source_sha256)
        for digest in self.prediction_lineage_sha256.values(): require_sha(digest)
        if set(self.assessment_components) != set(self.assessment_heads):
            raise ValueError('Identical row masks required for head/component admission')
        component_heads={}
        for rid,component in self.assessment_components.items():
            head=self.assessment_heads[rid]
            if not all(isinstance(v,str) and v for v in (rid,component,head)):
                raise ValueError('Nonempty row/component/head identities required')
            if component in component_heads and component_heads[component]!=head:
                raise ValueError('One dependency component cannot span original held-out heads')
            component_heads[component]=head
        for head,components in self.fitting_components_by_head.items():
            if (not isinstance(components,(tuple,list)) or len(set(components))!=len(components)
                    or any(not isinstance(v,str) or not v for v in components)):
                raise ValueError('Explicit unique fitting dependency components required')
        if any(type(v) is not int or v not in (-1,1) for v in self.prediction_directions.values()):
            raise ValueError('Prediction direction must explicitly map to native readout direction')
        for name in ('assessment_components','assessment_heads','fitting_components_by_head',
                     'prediction_lineage_sha256','arm_interpretations','prediction_directions'):
            object.__setattr__(self,name,_frozen(getattr(self,name)))


@dataclass(frozen=True)
class SubstitutionPrediction:
    record_id: str
    assay_id: str
    construct_id: str
    reference_sha256: str
    heldout_head: str
    component_id: str
    site_1based: int
    ref_aa: str
    alt_aa: str
    scores: Mapping[str,float | None]
    prediction_origin: str = 'frozen_runtime'

    def __post_init__(self):
        require_sha(self.reference_sha256)
        if type(self.site_1based) is not int or self.site_1based <= 0:
            raise ValueError('Positive reference residue numbering required')
        if (len(self.ref_aa)!=1 or len(self.alt_aa)!=1 or self.ref_aa not in 'ACDEFGHIKLMNPQRSTVWY'
                or self.alt_aa not in 'ACDEFGHIKLMNPQRSTVWY' or self.ref_aa == self.alt_aa):
            raise ValueError('Directed nonidentity amino acid substitution required')
        object.__setattr__(self,'scores',_scores(self.scores))


@dataclass(frozen=True)
class EndpointObservation:
    record_id: str
    value: float | None
    assay_id: str
    construct_id: str
    reference_sha256: str
    heldout_head: str
    component_id: str
    source_sha256: str
    native_units: str
    uncertainty: float | None = None

    def __post_init__(self):
        require_sha(self.reference_sha256)
        require_sha(self.source_sha256)
        if not self.native_units: raise ValueError('Explicit native measurement units required')
        if self.uncertainty is not None and (not _number(self.uncertainty) or self.uncertainty<0):
            raise ValueError('Native measurement uncertainty must be finite nonnegative or missing')
        if self.value is not None and not _number(self.value): raise ValueError('Native finite measurement or explicit missing required')


@dataclass(frozen=True)
class EndpointContract:
    endpoint_id: str
    source_sha256: str
    assay_id: str
    construct_id: str
    reference_sha256: str
    native_units: str
    target_definition: str
    admission_state: str
    release_ref: str
    observations: tuple[EndpointObservation,...]
    native_arm: str = 'native'
    uncertainty_definition: str = 'not_quantified'

    def __post_init__(self):
        require_sha(self.source_sha256); require_sha(self.reference_sha256)
        if not all((self.endpoint_id,self.assay_id,self.construct_id,self.native_units,self.target_definition,self.release_ref)):
            raise ValueError('Source/construct/native endpoint and release declarations required')
        object.__setattr__(self,'observations',tuple(self.observations))


@dataclass(frozen=True)
class PairSummary:
    total_pairs: int
    supported_pairs: int
    missing_pairs: int
    prediction_ties: int
    measured_ties: int
    directional_pairs: int
    concordance: float | None


@dataclass(frozen=True)
class BiologicalEvaluation:
    status: str
    reason: str
    endpoint_id: str
    native_units: str
    target_definition: str
    test_A: Mapping = field(default_factory=dict)
    test_B: Mapping = field(default_factory=dict)
    native_common: Mapping = field(default_factory=dict)
    singleton_groups: Mapping = field(default_factory=dict)
    group_support: tuple = ()
    arm_interpretations: Mapping = field(default_factory=dict)
    qualification: str = 'engineering_accounting_only'


def _aligned(predictions, contract):
    ids = [row.record_id for row in predictions]
    if len(set(ids)) != len(ids) or set(ids) != set(contract.assessment_components):
        raise ValueError('Predictions must exactly cover the admitted row mask without duplicates')
    arms = set(contract.prediction_lineage_sha256)
    for row in predictions:
        if (not row.record_id or not row.component_id or not row.heldout_head
                or row.component_id != contract.assessment_components[row.record_id]
                or row.heldout_head != contract.assessment_heads[row.record_id]):
            raise ValueError('Row component or original held-out head mismatch')
        if row.heldout_head not in contract.fitting_components_by_head:
            raise ValueError('Missing fitting exposure declaration for held-out head')
        # Each head must exclude all its assessment components. Components
        # assessed by a different original head may belong to its fitting set.
        heldout_components = {contract.assessment_components[rid] for rid,head in contract.assessment_heads.items()
                              if head == row.heldout_head}
        if heldout_components.intersection(contract.fitting_components_by_head[row.heldout_head]):
            raise ValueError('Training/assessment dependency component leakage')
        if set(row.scores) != arms: raise ValueError('Every declared prediction arm must retain missing rows explicitly')
        if row.prediction_origin not in {'frozen_runtime','component_cross_fitted_teacher'}:
            raise ValueError('Observed effects or in-sample teacher predictions cannot enter evaluation predictions')
    if set(contract.arm_interpretations) - arms: raise ValueError('Unbound arm interpretation')
    return tuple(sorted(arms))


def _pair_summary(pairs, arm, observations, *, native_arm=None, direction=1):
    total = len(pairs); supported = ties = measured_ties = directional = 0; credit = 0.
    for left,right in pairs:
        y1,y2 = observations[left.record_id].value,observations[right.record_id].value
        p1,p2 = left.scores[arm],right.scores[arm]
        if y1 is None or y2 is None or p1 is None or p2 is None: continue
        if native_arm is not None and (left.scores.get(native_arm) is None or right.scores.get(native_arm) is None): continue
        supported += 1
        tied = abs(p1-p2) <= 1e-12
        ties += int(tied)
        if y1 == y2:
            measured_ties += 1
            continue
        directional += 1
        credit += .5 if tied else float(((p1-p2)*direction>0) == (y1>y2))
    return PairSummary(total,supported,total-supported,ties,measured_ties,directional,
                       credit/directional if directional else None)


def evaluate_substitution_context(predictions, admitted_endpoint, split_contract):
    """TEST A: distinct ALTs at a site; TEST B: same REF→ALT across sites.

    Pairs never cross assay, construct, reference or original held-out heads.
    A measured tie has no directional denominator; a prediction tie gets 1/2
    credit on a measured unequal pair. Both remain explicitly counted.
    """
    rows = tuple(predictions); endpoint = admitted_endpoint; contract = split_contract
    if endpoint.admission_state != 'admitted':
        if any(obs.value is not None for obs in endpoint.observations):
            raise ValueError('Numerical targets lack admitted endpoint release')
        return BiologicalEvaluation('source_unavailable','Endpoint is not admitted for numerical evaluation',
                                    endpoint.endpoint_id,endpoint.native_units,endpoint.target_definition)
    arms = _aligned(rows,contract)
    if set(contract.prediction_directions)!=set(arms):
        raise ValueError('Every biological prediction arm needs a source-bound native-readout direction')
    if endpoint.source_sha256 != contract.assessment_source_sha256:
        raise ValueError('Endpoint source differs from frozen source admission')
    observations = {obs.record_id:obs for obs in endpoint.observations}
    if len(observations) != len(endpoint.observations) or set(observations) != {r.record_id for r in rows}:
        raise ValueError('Endpoint observations must exactly align with the admitted row IDs')
    for row in rows:
        obs = observations[row.record_id]
        if (row.assay_id,row.construct_id,row.reference_sha256) != (endpoint.assay_id,endpoint.construct_id,endpoint.reference_sha256):
            raise ValueError('Prediction source/construct/reference differs from qualified endpoint')
        if (obs.assay_id,obs.construct_id,obs.reference_sha256,obs.heldout_head,obs.component_id) != (
            row.assay_id,row.construct_id,row.reference_sha256,row.heldout_head,row.component_id):
            raise ValueError('Endpoint source/construct/head/component mismatch')
        if obs.source_sha256 != endpoint.source_sha256 or obs.native_units != endpoint.native_units:
            raise ValueError('Observation source pin or native units mismatch')
    results = {}; common = {}; singletons = {}; support = []
    for test in ('A','B'):
        groups = {}
        for row in rows:
            scope = (row.assay_id,row.construct_id,row.reference_sha256,row.heldout_head)
            key = scope + ((row.site_1based,row.ref_aa) if test == 'A' else (row.ref_aa,row.alt_aa))
            groups.setdefault(key,[]).append(row)
        pairs = []
        singletons[test] = sum(len(group)==1 for group in groups.values())
        for key,group in sorted(groups.items()):
            selected = [(a,b) for a,b in combinations(group,2)
                        if (a.alt_aa != b.alt_aa if test == 'A' else a.site_1based != b.site_1based)]
            pairs.extend(selected)
            support.append({'test':test,'group':key,'rows':len(group),'pairs':len(selected),
                            'missing_measurements':sum(observations[r.record_id].value is None for r in group)})
        results[test] = {arm:_pair_summary(pairs,arm,observations,direction=contract.prediction_directions[arm]) for arm in arms}
        common[test] = {}
        if endpoint.native_arm in arms:
            for arm in arms:
                matched = [(a,b) for a,b in pairs if a.scores[arm] is not None and b.scores[arm] is not None
                           and a.scores[endpoint.native_arm] is not None and b.scores[endpoint.native_arm] is not None]
                common[test][arm] = {'candidate':_pair_summary(matched,arm,observations,direction=contract.prediction_directions[arm]),
                                    'native':_pair_summary(matched,endpoint.native_arm,observations,
                                                         direction=contract.prediction_directions[endpoint.native_arm])}
    return BiologicalEvaluation('evaluated','Admitted endpoint-specific frozen predictions only',endpoint.endpoint_id,
        endpoint.native_units,endpoint.target_definition,results['A'],results['B'],common,singletons,
        tuple(support),contract.arm_interpretations)


@dataclass(frozen=True)
class ClinicalObservation:
    record_id: str
    label: int
    component_id: str
    heldout_head: str
    groups: tuple[str,...] = ()

    def __post_init__(self):
        if type(self.label) is not int or self.label not in (0,1): raise ValueError('Binary admitted clinical label required')
        object.__setattr__(self,'groups',tuple(self.groups))


@dataclass(frozen=True)
class ClinicalPrediction:
    record_id: str
    heldout_head: str
    component_id: str
    scores: Mapping[str,float | None]
    prediction_origin: str = 'frozen_runtime'

    def __post_init__(self): object.__setattr__(self,'scores',_scores(self.scores,probability=True))


@dataclass(frozen=True)
class ClinicalCohort:
    cohort_id: str
    source_sha256: str
    role_map_sha256: str
    release_ref: str
    admission_state: str
    rows: tuple[ClinicalObservation,...]
    primary_group: str = 'conservative'
    confirmation: bool = False
    exposed_development: bool = True
    confirmation_lock_sha256: str | None = None
    unexposed_release_ref: str | None = None
    prespecified_groups: tuple[str,...] = ('conservative',)

    def __post_init__(self):
        require_sha(self.source_sha256); require_sha(self.role_map_sha256)
        if self.confirmation_lock_sha256 is not None: require_sha(self.confirmation_lock_sha256)
        if not self.cohort_id or not self.release_ref or not self.primary_group:
            raise ValueError('Explicit admitted cohort, release and primary-group identity required')
        object.__setattr__(self,'rows',tuple(self.rows))
        object.__setattr__(self,'prespecified_groups',tuple(self.prespecified_groups))
        if self.primary_group not in self.prespecified_groups or len(set(self.prespecified_groups))!=len(self.prespecified_groups):
            raise ValueError('Unique prespecified groups must contain the primary group')


@dataclass(frozen=True)
class OperatingPointContract:
    thresholds: Mapping[str,float]
    frozen_receipt_sha256: str
    selection_role: str = 'frozen_protocol'
    locked_before_assessment: bool = True
    evee_dev_threshold_receipt_sha256: str | None = None

    def __post_init__(self):
        require_sha(self.frozen_receipt_sha256)
        if self.evee_dev_threshold_receipt_sha256 is not None: require_sha(self.evee_dev_threshold_receipt_sha256)
        if (self.selection_role not in {'frozen_protocol','fitting'} or not self.locked_before_assessment
                or any(not _number(v) for v in self.thresholds.values())):
            raise ValueError('Finite frozen fitting/protocol operating points required; assessment selection prohibited')
        object.__setattr__(self,'thresholds',_frozen(self.thresholds))


@dataclass(frozen=True)
class ClinicalMetrics:
    total: int
    scored: int
    unscored: int
    positives: int
    negatives: int
    coverage: float | None
    true_positives: int | None
    false_positives: int | None
    true_negatives: int | None
    false_negatives: int | None
    errors: int | None
    sensitivity: float | None
    specificity: float | None
    balanced_accuracy: float | None
    auroc: float | None
    auprc: float | None
    support: str


@dataclass(frozen=True)
class IntervalEstimate:
    lower: float | None
    upper: float | None
    requested_draws: int = 5000
    supported_draws: int = 0
    seed: int = 1729
    interpretation: str = 'paired component percentile interval conditional on exposed models/data'


@dataclass(frozen=True)
class GateDecision:
    passed: bool
    reasons: tuple[str,...]


@dataclass(frozen=True)
class ClinicalEvaluation:
    scopes: Mapping
    transitions: Mapping
    intervals: Mapping
    development_gate: GateDecision
    confirmation_gate: GateDecision
    scope_intervals: Mapping = field(default_factory=dict)
    calibration: str = 'uncalibrated scores; no probability or deployment claim'
    qualification: str = 'engineering_accounting_only'


def _ranking(labels,scores):
    positive = sum(labels); negative = len(labels)-positive
    if not positive or not negative: return None,None
    # Group tied scores: AUC gives 1/2 credit, average precision includes the
    # entire threshold tie block. No arbitrary ordering within a tied block.
    grouped = {}
    for y,score in zip(labels,scores):
        counts = grouped.setdefault(score,[0,0]); counts[y] += 1
    negatives_below = 0; wins = 0.
    for score,(neg,pos) in sorted(grouped.items()):
        wins += pos*(negatives_below+.5*neg)
        negatives_below += neg
    tp = fp = 0; ap = 0.
    for score,(neg,pos) in sorted(grouped.items(),reverse=True):
        tp += pos; fp += neg
        ap += (pos/positive)*(tp/(tp+fp))
    return wins/(positive*negative),ap


def _clinical_metrics(rows,predictions,arm,cutoff):
    total = len(rows); positives = sum(r.label for r in rows); negatives = total-positives
    scored_rows = [r for r in rows if predictions[r.record_id].scores[arm] is not None]
    n = len(scored_rows)
    if not n:
        return ClinicalMetrics(total,0,total,positives,negatives,0. if total else None,
                               None,None,None,None,None,None,None,None,None,None,'empty_or_unscored')
    tp=fp=tn=fn=0
    labels=[]; scores=[]
    for row in scored_rows:
        score=predictions[row.record_id].scores[arm]; call=score>=cutoff
        labels.append(row.label); scores.append(score)
        tp += int(call and row.label==1); fp += int(call and row.label==0)
        tn += int(not call and row.label==0); fn += int(not call and row.label==1)
    sens = tp/positives if positives and tp+fn==positives else None
    spec = tn/negatives if negatives and fp+tn==negatives else None
    ba = (sens+spec)/2 if sens is not None and spec is not None else None
    auc,ap = _ranking(labels,scores) if n==total else (None,None)
    support = 'supported' if ba is not None else ('incomplete_coverage' if n!=total else 'one_class')
    return ClinicalMetrics(total,n,total-n,positives,negatives,n/total,tp,fp,tn,fn,fp+fn,sens,spec,ba,auc,ap,support)


def _transitions(rows,predictions,points,candidate,comparator):
    result = dict(corrected_fp=0,new_fp=0,corrected_fn=0,new_fn=0,unchanged=0,not_comparable=0,total=len(rows))
    for row in rows:
        left=predictions[row.record_id].scores[comparator]; right=predictions[row.record_id].scores[candidate]
        if left is None or right is None:
            result['not_comparable'] += 1; continue
        before = left>=points[comparator]; after = right>=points[candidate]
        if before==after: result['unchanged'] += 1
        elif row.label==0: result['corrected_fp' if before else 'new_fp'] += 1
        else: result['new_fn' if before else 'corrected_fn'] += 1
    return result


def _intervals(rows,predictions,points,candidate,comparators,*,component_universe=None):
    keys = ('ba_gain','fpr_increase','fnr_increase')
    unsupported = lambda: {key:IntervalEstimate(None,None) for key in keys}
    if not rows: return {arm:unsupported() for arm in comparators}
    components = sorted(component_universe if component_universe is not None else {row.component_id for row in rows})
    indices={c:i for i,c in enumerate(components)}
    arms = (candidate,*comparators)
    matrices = np.zeros((len(components),len(arms),4),dtype=np.int64)
    complete = {arm:True for arm in arms}
    for row in rows:
        for a,arm in enumerate(arms):
            score=predictions[row.record_id].scores[arm]
            if score is None: complete[arm]=False; continue
            call=score>=points[arm]
            cell = 0 if row.label==1 and call else 1 if row.label==1 else 2 if call else 3
            matrices[indices[row.component_id],a,cell] += 1
    # One set of component draws is shared by all arms and all paired metrics.
    rng=np.random.default_rng(1729)
    sampled=np.empty((5000,len(arms),4),dtype=np.int64)
    for draw in range(5000):
        selected=rng.integers(0,len(components),len(components))
        sampled[draw]=matrices[selected].sum(axis=0)
    result={}
    for a,arm in enumerate(comparators,1):
        if not complete[candidate] or not complete[arm]: result[arm]=unsupported(); continue
        base=sampled[:,a]; cand=sampled[:,0]
        valid=(base[:,0]+base[:,1]>0)&(base[:,2]+base[:,3]>0)&(cand[:,0]+cand[:,1]>0)&(cand[:,2]+cand[:,3]>0)
        n=int(valid.sum())
        if not n: result[arm]=unsupported(); continue
        base=base[valid].astype(float); cand=cand[valid].astype(float)
        fnr_base=base[:,1]/(base[:,0]+base[:,1]); fnr_cand=cand[:,1]/(cand[:,0]+cand[:,1])
        fpr_base=base[:,2]/(base[:,2]+base[:,3]); fpr_cand=cand[:,2]/(cand[:,2]+cand[:,3])
        series=(-(fnr_cand+fpr_cand-fnr_base-fpr_base)/2,fpr_cand-fpr_base,fnr_cand-fnr_base)
        result[arm]={key:IntervalEstimate(float(np.quantile(values,.025)),float(np.quantile(values,.975)),
                                        supported_draws=n) for key,values in zip(keys,series)}
    return result


def select_fitting_threshold(scores, labels, max_false_positives, *, role, selection_receipt_sha256):
    """Lowest fitting threshold satisfying an FP budget; ties stay indivisible."""
    require_sha(selection_receipt_sha256)
    if role!='fitting': raise ValueError('Only fitting rows may select an operating threshold')
    scores=tuple(scores); labels=tuple(labels)
    if (not scores or len(scores)!=len(labels) or any(not _number(v) or not 0<=v<=1 for v in scores)
            or any(type(v) is not int or v not in (0,1) for v in labels)
            or type(max_false_positives) is not int or max_false_positives<0 or not any(v==0 for v in labels)):
        raise ValueError('Complete fitting scores, negative support and a nonnegative FP budget required')
    for threshold in sorted(set(scores)):
        if sum(score>=threshold and label==0 for score,label in zip(scores,labels))<=max_false_positives:
            return threshold
    return math.nextafter(max(scores),math.inf)


def evaluate_clinical_increment(predictions, admitted_cohort, operating_points, component_contract,
                               *, candidate_arm='candidate',baseline_arm='baseline',evee_arm='evee',
                               availability_arm='availability',reassignment_arm='reassigned'):
    """Account for every admitted row and fixed operating point, without fitting.

    Full-cohort metrics are unsupported when coverage is incomplete. Bootstrap
    intervals resample dependency components 5,000 times with seed 1729 and do
    not erase exposed model selection. Gate decisions confer no registry state.
    """
    rows=tuple(predictions); cohort=admitted_cohort; contract=component_contract; points=operating_points.thresholds
    if cohort.admission_state!='admitted': raise ValueError('Clinical targets require exact admitted cohort release')
    arms=_aligned(rows,contract)
    if cohort.role_map_sha256!=contract.role_map_sha256: raise ValueError('Clinical role-map admission mismatch')
    if cohort.source_sha256!=contract.assessment_source_sha256:
        raise ValueError('Clinical source admission mismatch')
    if set(points)!=set(arms): raise ValueError('Frozen operating points must bind all exact prediction arms')
    required={candidate_arm,baseline_arm,evee_arm,availability_arm,reassignment_arm}
    if not required<=set(arms): raise ValueError('Candidate/baseline/EVEE/availability/reassignment arms required')
    outcomes={row.record_id:row for row in cohort.rows}
    if len(outcomes)!=len(cohort.rows) or set(outcomes)!={r.record_id for r in rows}:
        raise ValueError('Clinical outcomes must align exactly to admitted prediction IDs')
    for row in cohort.rows:
        if row.component_id!=contract.assessment_components[row.record_id] or row.heldout_head!=contract.assessment_heads[row.record_id]:
            raise ValueError('Clinical component/original-head mismatch')
        if not set(row.groups)<=set(cohort.prespecified_groups):
            raise ValueError('Clinical group was not prespecified in the admitted cohort contract')
    pred={r.record_id:r for r in rows}
    scopes={'full_cohort':cohort.rows,
            'primary_conservative':tuple(r for r in cohort.rows if cohort.primary_group in r.groups),
            'complement':tuple(r for r in cohort.rows if cohort.primary_group not in r.groups)}
    for group in sorted(set(cohort.prespecified_groups)-{cohort.primary_group}):
        scopes['group:'+group]=tuple(r for r in cohort.rows if group in r.groups)
        scopes['complement:'+group]=tuple(r for r in cohort.rows if group not in r.groups)
    metrics={scope:{arm:_clinical_metrics(group,pred,arm,points[arm]) for arm in arms} for scope,group in scopes.items()}
    transitions={scope:{arm:_transitions(group,pred,points,candidate_arm,arm) for arm in (baseline_arm,evee_arm)}
                 for scope,group in scopes.items()}
    comparators=tuple(dict.fromkeys((baseline_arm,evee_arm,availability_arm,reassignment_arm)))
    intervals=_intervals(cohort.rows,pred,points,candidate_arm,comparators)
    universe={row.component_id for row in cohort.rows}
    scope_intervals={'full_cohort':intervals}
    for scope,group in scopes.items():
        if scope!='full_cohort':
            scope_intervals[scope]=_intervals(group,pred,points,candidate_arm,comparators,component_universe=universe)
    full=metrics['full_cohort']; candidate=full[candidate_arm]; baseline=full[baseline_arm]
    dev=[]
    if len({points[arm] for arm in (candidate_arm,baseline_arm,availability_arm,reassignment_arm)})!=1:
        dev.append('Development advancement requires a common frozen candidate/baseline/control operating point')
    if any(full[arm].support!='supported' for arm in required): dev.append('Incomplete coverage or unsupported clinical class scope')
    else:
        if candidate.false_negatives>=baseline.false_negatives: dev.append('No reduction in baseline FN')
        if candidate.false_positives>baseline.false_positives: dev.append('New net FP increase')
        if candidate.balanced_accuracy<=baseline.balanced_accuracy: dev.append('No overall baseline BA improvement')
        if candidate.auroc<baseline.auroc: dev.append('AUROC regression')
        for arm in (availability_arm,reassignment_arm):
            ci=intervals[arm]['ba_gain']
            if candidate.balanced_accuracy<=full[arm].balanced_accuracy or ci.lower is None or ci.lower<=0:
                dev.append('No positive paired BA lower bound over '+arm)
    confirm=[]
    if not (cohort.confirmation and not cohort.exposed_development and cohort.confirmation_lock_sha256
            and cohort.unexposed_release_ref and operating_points.locked_before_assessment
            and operating_points.evee_dev_threshold_receipt_sha256):
        confirm.append('Qualified locked unexposed confirmation provenance required')
    comparator=full[evee_arm]
    if candidate.support!='supported' or comparator.support!='supported': confirm.append('Confirmation scope unsupported')
    else:
        ci=intervals[evee_arm]
        if candidate.balanced_accuracy-comparator.balanced_accuracy<.02: confirm.append('EVEE BA gain below 0.02')
        if ci['ba_gain'].lower is None or ci['ba_gain'].lower<=0: confirm.append('Paired BA lower bound is not positive')
        for name in ('fpr_increase','fnr_increase'):
            if ci[name].upper is None or ci[name].upper>.005: confirm.append(name+' upper bound exceeds 0.005')
    return ClinicalEvaluation(metrics,transitions,intervals,GateDecision(not dev,tuple(dev)),
                              GateDecision(not confirm,tuple(confirm)),scope_intervals)
