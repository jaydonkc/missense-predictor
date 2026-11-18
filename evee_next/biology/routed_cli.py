"""JSON CLI for the deterministic, target-free missense evidence pipeline."""
import argparse
from dataclasses import fields,replace
import json,os,sys,tempfile
from pathlib import Path
from evee_next.report import VariantRequest
from evee_next.reference import ReferenceResolver
from .admission import EvidenceRegistry,pinned_json
from .baseline import CachedBaselineProvider,CachedThenQualifiedBaseline,QualifiedCheckpointBaseline
from .evidence_store import EvidenceStore
from .execution import BoundedFetcher,ExecutionBudget,RunnerSpec,QualifiedRunner
from .fusion import FusionArtifact
from .routed import RoutedBiologyPipeline
from .routing import EvidenceDecision,EvidenceRequest,ExecutionPolicy,canonical_sha,report_dict


class UnavailableProvider:
    def __init__(self,provider_id,reason='provider_not_configured',status='unavailable'):
        self.provider_id,self.reason,self.status=provider_id,reason,status
    def plan(self,context):
        v=context.variant
        identity={'assembly':v.assembly if v else 'GRCh38'}
        if v:identity.update(chromosome=v.chromosome,position_1based=v.position_1based,ref=v.ref,alt=v.alt)
        return (EvidenceRequest(self.provider_id,canonical_sha({'provider':self.provider_id,'status':self.status,'reason':self.reason}),identity),)
    def obtain(self,request,store,policy):
        return EvidenceDecision(request,True,self.status,self.reason)


def _json(path):
    # Config/request bytes are supplied by the operator, with strict duplicate
    # rejection; referenced model/source payloads have separate explicit pins.
    data=Path(path).read_bytes()
    if len(data)>1024**2:raise ValueError('Configuration/request exceeds 1 MiB')
    import hashlib
    return pinned_json(path,hashlib.sha256(data).hexdigest())


def configure(config_path,*,mode='cache-only'):
    config_path=Path(config_path).resolve();config=_json(config_path);root=config_path.parent
    if not isinstance(config,dict) or config.get('schema')!='routed-missense-config/v1' or not set(config)<={'schema','baseline','reference','cache','providers','policy','fusion','registry'}:
        raise ValueError('Strict routed-missense-config/v1 required')
    def path(value):return (root/value).resolve()
    policy_data=dict(config.get('policy',{}));policy_data['mode']=mode
    policy=ExecutionPolicy(**policy_data);budget=ExecutionBudget(policy);fetcher=BoundedFetcher(policy,budget)
    store=EvidenceStore(path(config.get('cache','evidence-cache')),
        lock_timeout_seconds=min(30.,policy.task_seconds,policy.invocation_seconds))
    base=config.get('baseline');baseline=CachedBaselineProvider(path(base['bundle']),base['sha256']) if base else None
    if base:
        if not set(base)<={'bundle','sha256','feature_outputs'}:raise ValueError('Unknown baseline configuration fields')
        if base.get('feature_outputs'):
            from .retained_features import RetainedFeatureProvider
            f=base['feature_outputs']
            head=baseline._provider(f['checkpoint_fold']).head
            features=RetainedFeatureProvider(path(f['manifest']),f['sha256'],checkpoint=head)
            baseline=CachedThenQualifiedBaseline(baseline,QualifiedCheckpointBaseline(head,features))
    ref=config.get('reference')
    if ref:
        if set(ref)!={'bundle','genome_sources','fasta_directory','concordance_receipt_sha256','manifest_sha256','translation_sha256'}:
            raise ValueError('Reference configuration requires exact source paths/pins')
        resolver=ReferenceResolver(path(ref['bundle']),path(ref['genome_sources']),path(ref['fasta_directory']),
            concordance_receipt_sha256=ref['concordance_receipt_sha256'],manifest_sha256=ref['manifest_sha256'],translation_sha256=ref['translation_sha256'])
    elif baseline:resolver=baseline
    else:raise ValueError('A qualified reference resolver or archived baseline is required')
    providers=[];source_config=config.get('providers',{})
    if not isinstance(source_config,dict) or not set(source_config)<={'structure','annotations','rna','reference_cache','disulfide','named_effects'}:
        raise ValueError('Unknown evidence provider configuration')
    for name in ('annotations','rna','structure'):
        p=source_config.get(name)
        if p is None:
            providers.append(UnavailableProvider('annotation' if name=='annotations' else name));continue
        try:
            if name=='structure':
                if not set(p)<={'manifest','sha256','artifact_root','source_root'}:raise ValueError('Unknown structure configuration fields')
                from .structure_provider import StructureEvidenceProvider
                providers.append(StructureEvidenceProvider(path(p['manifest']),manifest_sha256=p['sha256'],
                    artifact_root=path(p['artifact_root']) if p.get('artifact_root') else None,
                    source_root=path(p['source_root']) if p.get('source_root') else None,fetcher=fetcher))
            elif name=='annotations':
                if not set(p)<={'manifest','sha256','artifact_root'}:raise ValueError('Unknown annotation configuration fields')
                from .annotation_provider import FunctionalAnnotationProvider
                providers.append(FunctionalAnnotationProvider(path(p['manifest']),manifest_sha256=p['sha256'],
                    artifact_root=path(p['artifact_root']) if p.get('artifact_root') else None,fetcher=fetcher))
            else:
                from .rna_provider import SpliceEvidenceProvider,SpliceRecipe,LocalSpliceRunnerAdapter
                from .spliceai_source import SpliceSourceManifest,LocalSpliceSource,IndexedSpliceSource
                manifest=SpliceSourceManifest.from_bytes(path(p['manifest']).read_bytes(),expected_sha256=p['sha256'])
                if not set(p)<={'manifest','sha256','scores','engine','runner','runner_sha256'}:raise ValueError('Unknown RNA configuration fields')
                recipe=SpliceRecipe.from_manifest(manifest)
                if p.get('engine','precomputed')=='local':
                    if not p.get('runner'):
                        providers.append(UnavailableProvider('rna.spliceai','qualified_local_spliceai_runner_or_weights_or_annotation_absent'));continue
                    spec=RunnerSpec(**p['runner'])
                    runner=LocalSpliceRunnerAdapter(QualifiedRunner(spec,policy),expected_sha256=p['runner_sha256'])
                    recipe=replace(recipe,engine='local',runner_id=runner.runner_id,runner_sha256=runner.runner_sha256)
                    providers.append(SpliceEvidenceProvider(recipe,runner=runner,budget=budget))
                elif p.get('engine','precomputed')=='precomputed':
                    source=LocalSpliceSource(manifest,path(p['scores'])) if p.get('scores') else IndexedSpliceSource(manifest,fetcher)
                    providers.append(SpliceEvidenceProvider(recipe,source=source,budget=budget))
                else:raise ValueError('Unknown RNA execution engine')
        except Exception as exc:
            providers.append(UnavailableProvider('annotation' if name=='annotations' else name,'provider_configuration_failed: '+str(exc),'failed'))
    if 'named_effects' in source_config:
        from .named_effect_provider import EFFECTS, NamedEffectProvider
        p = source_config['named_effects']
        for endpoint in EFFECTS:
            try:
                if set(p) != {'bundle', 'sha256'}:
                    raise ValueError('Named effects require a target-free bundle and manifest hash')
                providers.append(NamedEffectProvider(path(p['bundle']), p['sha256'], endpoint))
            except Exception as exc:
                providers.append(UnavailableProvider('effect.'+endpoint, 'provider_configuration_failed: '+str(exc),'failed'))
    if 'disulfide' in source_config:
        p = source_config['disulfide']
        try:
            from .disulfide_provider import DisulfideLossProvider
            if not set(p) <= {'manifest', 'sha256', 'artifact_root'}:
                raise ValueError('Unknown disulfide configuration fields')
            providers.append(DisulfideLossProvider(path(p['manifest']), manifest_sha256=p['sha256'],
                artifact_root=path(p['artifact_root']) if p.get('artifact_root') else None, fetcher=fetcher))
        except Exception as exc:
            providers.append(UnavailableProvider('annotation.disulfide','provider_configuration_failed: '+str(exc),'failed'))
    if 'reference_cache' in source_config:
        p=source_config['reference_cache']
        try:
            from .reference_cache_provider import CachedReferenceContextProvider
            if set(p)!={'bundle','sha256'}:raise ValueError('Reference context cache requires only exact bundle/pin')
            providers.append(CachedReferenceContextProvider(path(p['bundle']),p['sha256']))
        except Exception as exc:
            providers.append(UnavailableProvider('structure.cached','provider_configuration_failed: '+str(exc),'failed'))
    fusion=None
    if config.get('fusion'):
        if not config.get('registry'):raise ValueError('Fusion requires an independently pinned registry')
        r=config['registry'];f=config['fusion']
        fusion=FusionArtifact(path(f['path']),f['sha256'],EvidenceRegistry(path(r['path']),r['sha256']))
    return RoutedBiologyPipeline(resolver,baseline,providers,store=store,fusion=fusion),policy


def write_report(path,report):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    payload=(json.dumps(report_dict(report),sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'-',dir=path.parent)
    try:
        with os.fdopen(fd,'wb') as stream:
            stream.write(payload);stream.flush();os.fsync(stream.fileno())
        os.link(tmp,path)  # create-new publication: never overwrite preserved runs
    finally:os.unlink(tmp)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('request',type=Path);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--policy',choices=('cache-only','bounded-retrieval'),default='cache-only')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args(argv);pipeline=None
    try:
        raw=_json(args.request)
        if not isinstance(raw,dict) or not set(raw)<={f.name for f in fields(VariantRequest)}:
            raise ValueError('Only the variant request fields are accepted')
        request=VariantRequest(**raw);pipeline,policy=configure(args.config,mode=args.policy)
        report=pipeline.analyze(request,policy=policy)
        if args.output:write_report(args.output,report)
        else:print(json.dumps(report_dict(report),sort_keys=True,indent=2,allow_nan=False))
        return 0
    except (ValueError,OSError,TypeError,KeyError) as exc:
        print(json.dumps({'status':'configuration_or_request_failed','reason':str(exc)}),file=sys.stderr)
        return 2
    finally:
        if pipeline and hasattr(pipeline.resolver,'close'):pipeline.resolver.close()


if __name__=='__main__':raise SystemExit(main())
