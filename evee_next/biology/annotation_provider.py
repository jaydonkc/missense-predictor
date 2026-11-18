"""Exact-reference functional points from pinned native UniProt JSON.

The source manifest must qualify completeness before absence is asserted. These
annotations describe reference sites and do not predict substitution damage.
"""
from pathlib import Path
from urllib.parse import urlsplit

from .execution import BudgetExceeded,require_current_policy
from .routing import EvidenceDecision, EvidenceRequest, canonical_sha
from .structure_provider import _digest, _manifest, _path, _read, _request, _variant_budget_key, MAX_PAYLOAD

KINDS = {'Active site':'active_site_overlap', 'Binding site':'binding_site_overlap',
         'Site':'functional_site_overlap'}


class FunctionalAnnotationProvider:
    provider_id = 'annotation'

    def __init__(self, manifest_path, *, manifest_sha256, artifact_root=None, fetcher=None):
        self.manifest = _manifest(manifest_path, manifest_sha256, 'functional-annotation-sources/v1')
        self.manifest_sha256 = manifest_sha256
        self.artifact_root = Path(artifact_root or Path(manifest_path).parent).resolve()
        self.fetcher = fetcher
        self._raw = {}
        self.records = self.manifest.get('records')
        if not isinstance(self.records,tuple): raise ValueError('Qualified native annotation records are required')
        self.recipe_sha256 = canonical_sha({'provider':self.provider_id,'version':1,
            'manifest_sha256':manifest_sha256,'code_sha256':_digest(Path(__file__).read_bytes()),
            'shared_code_sha256':_digest(Path(__file__).with_name('structure_provider.py').read_bytes())})

    def plan(self, context):
        return (_request(self.provider_id,self.recipe_sha256,context),)

    def _failure(self,request,status,reason):
        return EvidenceDecision(request,True,status,reason,
            metadata={'annotation_state':'missing','reference_context_only':True})

    def _payload(self,record,request,store,policy):
        if store is None:return self._acquire_payload(record,request,None,policy)
        identity={k:request.identity[k] for k in ('protein_id','protein_sha256')}
        raw_request=EvidenceRequest('annotation.artifact',self.recipe_sha256,identity,
            {'sha256':record['sha256'],'accession':record['accession'],'entry_version':record['entry_version']})
        observed={}
        def producer():
            raw,disposition,execution=self._acquire_payload(record,request,store,policy)
            observed.update(disposition=disposition,execution=execution)
            return raw,EvidenceDecision(raw_request,True,'available','Hash-qualified native annotation bytes',
                {'payload_bytes':len(raw)},{'payload_bytes':'byte'},
                {'manifest_sha256':self.manifest_sha256,'artifact_sha256':record['sha256'],
                    'source_version':self.manifest['source_version'],'license':self.manifest['license']},
                {'artifact_only':True},qualification='artifact_integrity',execution=execution)
        store.get_or_compute(raw_request,producer)
        raw=store.read_artifact(raw_request)
        if raw is None:raise ValueError('Published shared annotation artifact is missing')
        self._raw[record['sha256']]=raw
        return raw,observed.get('disposition','raw_cache_hit'),observed.get('execution','not_run')

    def _acquire_payload(self, record, request, store, policy):
        pin = record['sha256']
        if pin in self._raw:return self._raw[pin], 'raw_memory_hit', 'not_run'
        if store is not None:
            data = store.get_artifact(pin)
            if data is not None:
                self._raw[pin] = data
                return data, 'raw_cache_hit', 'not_run'
        path = _path(self.artifact_root,record['path'])
        if path.exists():
            data = _read(path,pin)
            self._raw[pin] = data
            return data, 'local_payload', 'not_run'
        if policy.mode == 'cache-only': raise FileNotFoundError('Qualified annotation payload is absent')
        url = record.get('url')
        if self.fetcher is None or not isinstance(url,str):
            raise PermissionError('Annotation retrieval needs a configured bounded fetcher and source URL')
        parsed = urlsplit(url)
        if parsed.scheme != 'https' or parsed.hostname not in policy.allowed_hosts:
            raise PermissionError('Annotation host is not authorized by execution policy')
        require_current_policy(policy,fetcher=self.fetcher)
        data = self.fetcher.fetch(url,variant_key=_variant_budget_key(request),expected_etag=record.get('etag'))
        if len(data) > MAX_PAYLOAD or _digest(data) != pin:
            raise ValueError('Retrieved annotation payload size/SHA256 mismatch')
        self._raw[pin] = data
        return data, 'retrieved', 'bounded_retrieval'

    def _derive(self,request,store,policy):
        state = request.parameters.get('identity_status')
        if state != 'qualified':
            return b'', self._failure(request,'ambiguous' if state == 'ambiguous' else 'unavailable',
                'A unique validated missense mapping and exact full protein reference are required')
        records = [r for r in self.records if r.get('protein_id') == request.identity.get('protein_id')
            and r.get('protein_sha256') == request.identity.get('protein_sha256')]
        if not records: return b'', self._failure(request,'unavailable','No qualified exact-reference functional annotation')
        if len(records)>1:return b'',self._failure(request,'ambiguous','Multiple annotation sources need explicit selection')
        record = records[0]
        try:
            import json
            raw,disposition,execution = self._payload(record,request,store,policy)
            if _digest(raw) != record['sha256']: raise ValueError('Annotation payload SHA256 mismatch')
            data = json.loads(raw)
            if data['primaryAccession'] != record['accession'] or data['entryAudit']['entryVersion'] != record['entry_version']:
                raise ValueError('Native annotation accession/version differs from qualified source')
            sequence = data['sequence']['value']; p = request.identity['residue_1based']
            if _digest(sequence.encode()) != request.identity['protein_sha256'] or not 1 <= p <= len(sequence) or sequence[p-1] != request.identity['ref_aa']:
                return raw,self._failure(request,'unavailable','Annotation/reference full sequence correspondence failed')
            complete = record.get('complete_functional_points') is True
            overlaps = {name:0 if complete else None for name in KINDS.values()}
            uncertain = set()
            features = data.get('features')
            if not isinstance(features,list):
                return raw,self._failure(request,'unavailable','Native annotation feature collection is missing')
            for feature in features:
                if not isinstance(feature,dict):raise ValueError('Malformed native functional annotation record')
                name = KINDS.get(feature.get('type'))
                if name is None: continue
                location = feature.get('location',{});start = location.get('start',{});end = location.get('end',{})
                a,b = start.get('value'),end.get('value')
                exact = all(x.get('modifier') == 'EXACT' for x in (start,end))
                if not exact or type(a) is not int or type(b) is not int or not 1 <= a <= b <= len(sequence):
                    uncertain.add(name);continue
                # Only direct, exact single-residue points are admitted. A
                # functional region must not silently become point evidence.
                if a != b: uncertain.add(name);continue
                if p == a: overlaps[name] = 1
            for name in uncertain:
                if overlaps[name] != 1: overlaps[name] = None
            known = sum(v is not None for v in overlaps.values())
            if not known:return raw,self._failure(request,'unavailable','No qualified functional-point presence/absence at this site')
            annotation_state = 'overlap' if 1 in overlaps.values() else 'no_overlap' if known == len(overlaps) else 'incomplete'
            provenance = dict(manifest_sha256=self.manifest_sha256,artifact_sha256=record['sha256'],
                source_id=self.manifest['source_id'],source_version=self.manifest['source_version'],
                license=self.manifest['license'],accession=record['accession'],entry_version=str(record['entry_version']),
                sequence_source=self.manifest['sequence_reference']['source'],
                sequence_release=self.manifest['sequence_reference']['release'])
            return raw,EvidenceDecision(request,True,'available' if known == len(overlaps) else 'partial',
                'Exact reference functional-point annotations; no mutation-effect inference',overlaps,
                {name:'binary_direct_point_overlap' for name in overlaps},provenance,
                dict(annotation_state=annotation_state,reference_context_only=True,
                    supported_types=list(KINDS),complete_functional_points=complete,
                    uncertain_types=sorted(uncertain)),cache_disposition=disposition,execution=execution)
        except FileNotFoundError as exc:return b'',self._failure(request,'unavailable',str(exc))
        except PermissionError as exc:return b'',self._failure(request,'pending_authorization',str(exc))
        except BudgetExceeded as exc:return b'',self._failure(request,'unavailable',str(exc))
        except (KeyError,ValueError,TypeError,OSError,UnicodeError) as exc:
            return b'',self._failure(request,'failed',f'Qualified annotation extraction failed: {exc}')

    def obtain(self,request,store,policy):
        if request.provider_id != self.provider_id or request.recipe_sha256 != self.recipe_sha256:
            return self._failure(request,'failed','Request does not bind this provider source/recipe')
        if store is None:return self._derive(request,store,policy)[1]
        try: return store.get_or_compute(request,lambda:self._derive(request,store,policy))
        except (ValueError,OSError) as exc:
            return self._failure(request,'failed',f'Evidence cache integrity failure: {exc}')
