"""Inference-only frozen pre-score adapter. Never trains or admits itself."""
from dataclasses import replace
import math
from pathlib import Path
import numpy as np
from scipy.special import expit
from .admission import pinned_json
from .routing import require_sha,canonical_sha,validate_request_identity,_frozen
from evee_next.interpretation.neural import immutable


class FusionArtifact:
    def __init__(self,path,artifact_sha256,registry):
        self.path=Path(path);self.sha256=artifact_sha256
        self.artifact=pinned_json(path,artifact_sha256);a=self.artifact
        expected={'schema','checkpoint','feature_map_id','cutoff','features','R','qualification'}
        if set(a)!=expected or a['schema']!='routed-frozen-adapter/v1': raise ValueError('Frozen adapter schema required')
        require_sha(a['feature_map_id'])
        q=a['qualification']
        if set(q)!={'path','sha256'}:raise ValueError('Qualification receipt binding required')
        receipt=pinned_json(self.path.parent/q['path'],q['sha256'])
        if receipt.get('schema')!='routed-evidence-qualification/v1' or receipt.get('state')!='confirmed_scoring' or receipt.get('use')!='clinical_increment':
            raise ValueError('Confirmed clinical-use receipt required')
        self.entry=registry.require_confirmed(artifact_sha256,a['checkpoint'],a['feature_map_id'],q['sha256'])
        core={k:v for k,v in a.items() if k!='qualification'}
        scope={'assembly':self.entry['assembly'],'protein_ids':self.entry['protein_ids']}
        if (receipt.get('checkpoint'),receipt.get('feature_map_id'),receipt.get('adapter_payload_sha256'),receipt.get('scope')) != (
            a['checkpoint'],a['feature_map_id'],canonical_sha(core),scope):
            raise ValueError('Qualification receipt does not bind exact adapter/features/recipes/scope')
        evaluation=receipt.get('evaluation',{})
        if set(evaluation)!={'path','sha256'}:raise ValueError('Pinned independent evaluation receipt required')
        evaluated=pinned_json(self.path.parent/evaluation['path'],evaluation['sha256'])
        if (evaluated.get('schema'),evaluated.get('status'),evaluated.get('use'),evaluated.get('adapter_payload_sha256'),evaluated.get('scope')) != (
            'routed-clinical-confirmation/v1','independently_confirmed','clinical_increment',canonical_sha(core),scope):
            raise ValueError('Evaluation receipt does not qualify this exact clinical increment')
        if any(evaluated.get(k) is not True for k in ('model_locked','operating_points_locked','cohort_unexposed')):
            raise ValueError('Independent locked confirmation is required for routine scoring')
        if not a['features']:raise ValueError('Explicit ordered adapter features required')
        seen=set()
        for f in a['features']:
            if set(f)!={'provider_id','recipe_sha256','field','unit','kind','mean','scale'}:
                raise ValueError('Explicit feature/scaler contract required')
            require_sha(f['recipe_sha256'])
            identity=(f['provider_id'],f['recipe_sha256'],f['field'],f['kind'])
            if identity in seen:raise ValueError('Duplicate adapter feature')
            seen.add(identity)
            if f['kind'] not in {'value','availability'} or not f['unit']:
                raise ValueError('Unsupported scoring feature type')
            if any(isinstance(f[k],bool) or not isinstance(f[k],(int,float)) or not math.isfinite(f[k]) for k in ('mean','scale')) or f['scale']<=0:
                raise ValueError('Finite fitting-only scaler required')
            if f['kind']=='availability' and (f['mean'],f['scale'])!=(0,1):
                raise ValueError('Availability flags remain literal 0/1')
        self.R=np.asarray(a['R'],dtype=np.float64)
        if self.R.ndim!=2 or self.R.shape[0]!=len(a['features']) or not np.isfinite(self.R).all():
            raise ValueError('Finite feature-by-hidden adapter matrix required')
        self.R=immutable(self.R)
        self.artifact=_frozen(a);self.entry=_frozen(self.entry)

    def apply(self,baseline,decisions,context):
        a=self.artifact;state=baseline.state
        if state is None or baseline.prediction.score is None:raise ValueError('Complete base required')
        if context.variant is None or context.mapping is None or state.variant_key!=context.variant.variant_id:
            raise ValueError('Fusion state differs from variant identity')
        if (state.checkpoint,state.feature_map_id)!=(a['checkpoint'],a['feature_map_id']):
            raise ValueError('Fusion checkpoint/feature map differs')
        if self.entry['assembly']!=context.variant.assembly or (self.entry['protein_ids'] and context.mapping.protein_id not in self.entry['protein_ids']):
            raise ValueError('Variant outside scoring qualification scope')
        if a['cutoff']!=baseline.prediction.cutoff:raise ValueError('Adapter changes locked common operating point')
        if self.R.shape[1]!=len(state.a0):raise ValueError('Adapter hidden dimension differs')
        z=[];used={}
        for f in a['features']:
            matches=[(i,d) for i,d in enumerate(decisions) if d.request.provider_id==f['provider_id'] and d.request.recipe_sha256==f['recipe_sha256']]
            matches=[(i,d) for i,d in matches if d.request.identity.get('transcript_id') in {None,context.mapping.transcript_id}]
            if len(matches)!=1:return baseline.prediction,decisions
            i,d=matches[0]
            validate_request_identity(d.request,context,selected=True)
            if f['kind']=='availability':value=float(d.status in {'available','partial'} and d.values.get(f['field']) is not None)
            else:
                if d.status not in {'available','partial'} or d.values.get(f['field']) is None:return baseline.prediction,decisions
                if d.units.get(f['field'])!=f['unit']:raise ValueError('Scoring evidence unit differs')
                value=d.values[f['field']];used.setdefault(i,set()).add(f['field'])
            z.append((value-f['mean'])/f['scale'])
        # No value information means no evidence-based update. This also ensures
        # missingness alone cannot turn unavailable evidence into apparent damage.
        if not used:return baseline.prediction,decisions
        hidden=np.asarray(state.a0)+np.asarray(z)@self.R
        logit=state.b1+math.fsum(float(x*w) for x,w in zip(np.maximum(hidden,0),state.W1))
        if not math.isfinite(logit):raise ValueError('Nonfinite adapter output')
        final=replace(baseline.prediction,score=float(expit(logit)),scope='confirmed_scoring',reason='qualified pre-score evidence adapter '+self.sha256)
        consumed=tuple(replace(d,scoring_used=tuple(sorted(used.get(i,())))) for i,d in enumerate(decisions))
        return final,consumed
