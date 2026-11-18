"""Independent deterministic evidence routes around a standalone baseline."""
from .routing import (BaselinePrediction, BaselineResult, EvidenceDecision, EvidenceRequest,
    ExecutionPolicy, RoutedBiologyReport, VariantContext, canonical_sha,validate_request_identity)


class RoutedBiologyPipeline:
    def __init__(self, resolver, baseline=None, providers=(), *, store=None, fusion=None):
        self.resolver = resolver
        self.baseline = baseline
        self.providers = tuple(sorted(providers, key=lambda p: p.provider_id))
        if len({p.provider_id for p in self.providers}) != len(self.providers):
            raise ValueError('Unique providers required')
        self.store, self.fusion = store, fusion
        if fusion is not None:
            from .fusion import FusionArtifact
            if not isinstance(fusion,FusionArtifact):raise ValueError('A hash-bound admitted frozen fusion artifact is required')

    def analyze(self, request, *, policy=None):
        policy = policy or ExecutionPolicy()
        context = (self.resolver.context_for(request) if hasattr(self.resolver,'context_for')
                   else VariantContext.from_request(request,self.resolver))
        base = BaselineResult(BaselinePrediction(context.status, reason=context.detail))
        if context.genomic_validated and context.mapping is not None and self.baseline is not None:
            try: base = self.baseline.predict(context)
            except Exception as exc:
                base = BaselineResult(BaselinePrediction('failed',reason=str(exc)))
        if not isinstance(base,BaselineResult): raise TypeError('Typed baseline result required')
        decisions = []
        for provider in self.providers:
            try:
                requests = tuple(provider.plan(context))
                if not requests: raise ValueError('Provider must record why its route is skipped')
                if any(r.provider_id != provider.provider_id for r in requests):
                    raise ValueError('Provider request ID differs')
                if context.genomic_validated:
                    for item in requests:validate_request_identity(item,context)
            except Exception as exc:
                fallback = EvidenceRequest(provider.provider_id, canonical_sha({'provider':provider.provider_id,'planning_failure':True}),
                    {'assembly':context.variant.assembly if context.variant else 'GRCh38'})
                decisions.append(EvidenceDecision(fallback,False,'failed','planning_failed: '+str(exc)))
                continue
            for item in requests:
                if not context.genomic_validated:
                    result = EvidenceDecision(item,False,'not_applicable','unqualified_genomic_identity: '+context.status)
                else:
                    try:
                        result = provider.obtain(item,self.store,policy)
                        if not isinstance(result,EvidenceDecision) or result.request != item:
                            raise ValueError('Evidence result does not bind its exact request')
                        if result.scoring_used: raise ValueError('A provider cannot declare itself consumed in scoring')
                    except Exception as exc:
                        result = EvidenceDecision(item,True,'failed','provider_failed: '+str(exc))
                decisions.append(result)
        final, fusion_status = base.prediction, 'not_configured'
        if self.fusion is not None:
            try:
                if base.state is None: raise ValueError('Compatible complete baseline state unavailable')
                final, decisions = self.fusion.apply(base,tuple(decisions),context)
                fusion_status = 'applied' if final is not base.prediction else 'baseline_inherited'
            except Exception as exc:
                final, fusion_status = base.prediction, 'rejected: '+str(exc)
        return RoutedBiologyReport(context,base.prediction,final,tuple(decisions),fusion_status)


if __name__ == '__main__':
    from .routed_cli import main
    raise SystemExit(main())
