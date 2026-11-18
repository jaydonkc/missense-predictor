"""Source-pinned reference context; no mutation-effect or clinical inference.

Only the fixed AFDB v6 reviewed parser and C-alpha recipe below are executable.
The manifest is a caller-supplied trust anchor, not a discovered download receipt.
"""
import ast
from contextlib import nullcontext
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import zlib
from threading import RLock
from types import MappingProxyType
from urllib.parse import urlsplit

import numpy as np

from .execution import BudgetExceeded,require_current_policy
from .routing import EvidenceDecision, EvidenceRequest, canonical_sha, require_sha

PARSER_SHA256 = '964dc6b27d2509899ea59e652f526153bf4bba2d873327cca63445a9d7dcb8a1'
CIF_SHA256 = '0b418606f5fea90a4841f54b9415ebd7351677f0a670c55e9c4ec4aaf72a3d8a'
GEOMETRY_SHA256 = 'd27c7499c8176aa7e2fd2412a00efe16cc076052f1746f8c43823be257c1f937'
MAX_PAYLOAD = 64 * 1024**2
CANONICAL = frozenset('ACDEFGHIKLMNPQRSTVWY')
GEOMETRY_NAMES = ('CA_contacts_8A_sep_gt4', 'CA_contacts_12A_sep_gt4',
    'CA_contacts_12A_sep_gt12', 'nearest_nonlocal_CA_distance_A',
    'mean_5_nearest_nonlocal_CA_distance_A', 'CA_distance_from_chain_centroid_A')


def _immutable(value):
    if isinstance(value,dict): return MappingProxyType({k:_immutable(v) for k,v in value.items()})
    if isinstance(value,list): return tuple(_immutable(v) for v in value)
    return value


def _digest(data): return hashlib.sha256(data).hexdigest()


def _path(root, relative):
    value = Path(relative)
    if value.is_absolute() or '..' in value.parts:
        raise ValueError('Source paths must be contained relative paths')
    target = (root / value).resolve()
    if not target.is_relative_to(root.resolve()):
        raise ValueError('Source path escapes the qualified root')
    return target


def _read(path, expected, maximum=MAX_PAYLOAD):
    require_sha(expected)
    with path.open('rb') as stream:
        data = stream.read(maximum + 1)
    if len(data) > maximum: raise ValueError('Source exceeds bounded payload size')
    if _digest(data) != expected: raise ValueError('Source payload SHA256 mismatch')
    return data


def _manifest(path, pin, schema):
    data = json.loads(_read(Path(path), pin))
    if data.get('schema') != schema or data.get('qualification') != 'reference_context':
        raise ValueError('An explicitly qualified reference-context source manifest is required')
    for key in ('source_id', 'source_version', 'license'):
        if not isinstance(data.get(key), str) or not data[key].strip():
            raise ValueError('Versioned source identity and license are required')
    seq = data.get('sequence_reference', {})
    if not isinstance(seq, dict) or not all(isinstance(seq.get(k), str) and seq[k] for k in ('source','release')):
        raise ValueError('An explicit sequence reference source and release are required')
    return _immutable(data)


def _safe_definitions(path, pin, names, namespace):
    """Compile named, immutable recipe definitions without module entrypoints.

    Source digest is compared to a hardcoded reviewed implementation before any
    execution. This is not a facility for executing arbitrary manifest code.
    """
    tree = ast.parse(_read(path, pin), filename=str(path))
    chosen = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names:
            chosen.append(node)
        elif isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in names for t in node.targets):
            chosen.append(node)
    provided = {n.name if isinstance(n, (ast.FunctionDef, ast.ClassDef)) else n.targets[0].id for n in chosen}
    if provided != set(names): raise ValueError('Pinned recipe definitions differ')
    exec(compile(ast.Module(body=chosen, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def _variant_budget_key(request):
    return canonical_sha({key:request.identity[key] for key in
        ('assembly','chromosome','position_1based','ref','alt') if key in request.identity})


def _request(provider, recipe, context):
    v = context.variant
    identity = {} if v is None else dict(assembly=v.assembly, organism=v.organism,
        chromosome=v.chromosome, position_1based=v.position_1based, ref=v.ref, alt=v.alt)
    m = context.mapping
    if m is None:
        return EvidenceRequest(provider, recipe, identity or {'organism':'Homo sapiens'},
            {'identity_status':'ambiguous' if context.genomic_validated else 'unavailable'})
    identity.update(transcript_id=m.transcript_id, protein_id=m.protein_id, strand=m.strand,
        residue_1based=m.residue_1based, ref_aa=m.ref_aa, alt_aa=m.alt_aa)
    transcripts = [t for t in context.transcripts if t.transcript_id == m.transcript_id and t.protein_id == m.protein_id]
    sequence = getattr(transcripts[0], 'protein_sequence', None) if len(transcripts) == 1 else None
    status = 'qualified' if context.genomic_validated and isinstance(sequence, str) and sequence and set(sequence) <= CANONICAL else 'unavailable'
    if status == 'qualified':
        if m.residue_1based > len(sequence) or sequence[m.residue_1based-1] != m.ref_aa:
            status = 'unavailable'
        else: identity['protein_sha256'] = _digest(sequence.encode())
    return EvidenceRequest(provider, recipe, identity, {'identity_status':status})


class StructureEvidenceProvider:
    provider_id = 'structure'
    GEOMETRY_NAMES = GEOMETRY_NAMES

    def __init__(self, manifest_path, *, manifest_sha256, artifact_root=None, source_root=None, fetcher=None):
        self.manifest = _manifest(manifest_path, manifest_sha256, 'structure-evidence-sources/v1')
        self.manifest_sha256 = manifest_sha256
        self.artifact_root = Path(artifact_root or Path(manifest_path).parent).resolve()
        self.source_root = Path(source_root or Path(__file__).resolve().parents[2]).resolve()
        self.fetcher = fetcher
        self._raw = {}
        self._lock = RLock()
        self._geometry = None
        self._parser = None
        self._cif = None
        models = self.manifest.get('models')
        if not isinstance(models, tuple): raise ValueError('Qualified model records are required')
        self.models = models
        self.recipe_sha256 = canonical_sha({'provider': self.provider_id, 'version':1,
            'manifest_sha256':manifest_sha256, 'code_sha256':_digest(Path(__file__).read_bytes())})

    def plan(self, context):
        return (_request(self.provider_id, self.recipe_sha256, context),)

    def _source(self, kind, expected):
        entry = self.manifest.get('recipes', {}).get(kind, {})
        if entry.get('sha256') != expected:
            raise ValueError('No applicable pinned ' + kind + ' recipe')
        return _path(self.source_root, entry['path'])

    def geometry_recipe(self):
        with self._lock:
            if self._geometry is None:
                self._geometry = _safe_definitions(self._source('geometry', GEOMETRY_SHA256), GEOMETRY_SHA256,
                    {'NAMES','descriptors'}, {'np':np})
            return self._geometry

    def _parse(self, raw, model):
        with self._lock:
            if self._parser is None:
                self._parser = _safe_definitions(self._source('parser', PARSER_SHA256), PARSER_SHA256,
                    {'AA3','MEMBER_RE','single','clean_seq','rows','parse_record'},
                    {'Path':Path,'re':re,'np':np,'hashlib':hashlib,'json':json})
            if self._cif is None:
                # The pinned tokenizer needs only an already-open StringIO. No
                # vendored package path insertion or dependency install occurs.
                self._cif = _safe_definitions(self._source('cif_dictionary', CIF_SHA256), CIF_SHA256,
                    {'MMCIF2Dict'}, {'as_handle':nullcontext})['MMCIF2Dict']
        fmt = model.get('format')
        if fmt == 'mmcif.gz':
            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream: decoded = stream.read(MAX_PAYLOAD+1)
            if len(decoded) > MAX_PAYLOAD: raise ValueError('Bounded CIF decompression exceeded')
        elif fmt == 'mmcif': decoded = raw
        else: raise ValueError('Only qualified mmCIF payloads are supported')
        if _digest(decoded) != model.get('raw_sha256'): raise ValueError('Raw mmCIF SHA256 mismatch')
        dictionary = self._cif(io.StringIO(decoded.decode('utf-8')))
        return self._parser['parse_record'](dictionary, model['member'])

    def _provenance(self, model):
        return dict(manifest_sha256=self.manifest_sha256, artifact_sha256=model['sha256'],
            source_id=self.manifest['source_id'], source_version=self.manifest['source_version'],
            license=self.manifest['license'], parser_sha256=PARSER_SHA256,
            cif_dictionary_sha256=CIF_SHA256,
            sequence_source=self.manifest['sequence_reference']['source'],
            sequence_release=self.manifest['sequence_reference']['release'],model_id=model['model_id'])

    def _failure(self, request, status, reason, **kw):
        return EvidenceDecision(request, True, status, reason, metadata={'reference_context_only':True}, **kw)

    def _payload(self, model, request, store, policy):
        if store is None: return self._acquire_payload(model,request,None,policy)
        # A model-level key shares the raw retrieval lock across residue/ALT
        # requests and provider instances. Derived evidence keeps its own key.
        identity = {k:request.identity[k] for k in ('protein_id','protein_sha256')}
        identity.update(chain_id=model['chain_id'],fragment_id=model['model_id'])
        raw_request = EvidenceRequest('structure.artifact',self.recipe_sha256,identity,
            {'sha256':model['sha256'],'source_version':self.manifest['source_version']})
        observed = {}
        def producer():
            raw,disposition,execution = self._acquire_payload(model,request,store,policy)
            observed.update(disposition=disposition,execution=execution)
            return raw,EvidenceDecision(raw_request,True,'available','Hash-qualified reference artifact bytes',
                {'payload_bytes':len(raw)},{'payload_bytes':'byte'},self._provenance(model),
                {'artifact_only':True},qualification='artifact_integrity',execution=execution)
        store.get_or_compute(raw_request,producer)
        raw = store.read_artifact(raw_request)
        if raw is None: raise ValueError('Published shared artifact is missing')
        self._raw[model['sha256']] = raw
        return raw,observed.get('disposition','raw_cache_hit'),observed.get('execution','not_run')

    def _acquire_payload(self, model, request, store, policy):
        pin = model['sha256']; require_sha(pin)
        if pin in self._raw: return self._raw[pin], 'raw_memory_hit', 'not_run'
        if store is not None:
            try: cached = store.get_artifact(pin)
            except FileNotFoundError: cached = None
            if cached is not None:
                if _digest(cached) != pin: raise ValueError('Cached raw artifact SHA256 mismatch')
                self._raw[pin] = cached
                return cached, 'raw_cache_hit', 'not_run'
        path = _path(self.artifact_root, model['path'])
        if path.exists():
            raw = _read(path, pin)
            self._raw[pin] = raw
            return raw, 'local_payload', 'not_run'
        if policy.mode == 'cache-only': raise FileNotFoundError('Qualified local structure payload is absent; a historical receipt is not payload')
        url = model.get('url')
        if self.fetcher is None or not isinstance(url,str):
            raise PermissionError('Individual-model retrieval requires a configured bounded fetcher and source URL')
        parsed = urlsplit(url)
        if parsed.scheme != 'https' or parsed.hostname not in policy.allowed_hosts:
            raise PermissionError('Individual-model host is not authorized by execution policy')
        require_current_policy(policy,fetcher=self.fetcher)
        raw = self.fetcher.fetch(url, variant_key=_variant_budget_key(request), expected_etag=model.get('etag'))
        if len(raw) > MAX_PAYLOAD or _digest(raw) != pin: raise ValueError('Retrieved source payload size/SHA256 mismatch')
        self._raw[pin] = raw
        return raw, 'retrieved', 'bounded_retrieval'

    def _derive(self, request, store, policy):
        if request.parameters.get('identity_status') != 'qualified':
            state = request.parameters.get('identity_status')
            return b'', self._failure(request, 'ambiguous' if state == 'ambiguous' else 'unavailable',
                'A unique validated missense mapping and exact full protein reference are required')
        matches = [m for m in self.models if m.get('protein_id') == request.identity.get('protein_id')
            and m.get('protein_sha256') == request.identity.get('protein_sha256')]
        if not matches: return b'', self._failure(request, 'unavailable', 'No qualified exact-reference structure model')
        if len(matches) > 1: return b'', self._failure(request, 'ambiguous', 'Multiple qualified models require explicit source selection')
        model = matches[0]
        if model.get('source_type') != 'alphafold_db_v6':
            return b'', self._failure(request, 'unavailable', 'Alternate model source/confidence interpretation is unqualified')
        sequence = model.get('reference_sequence')
        if not isinstance(sequence,str) or _digest(sequence.encode()) != request.identity.get('protein_sha256'):
            return b'', self._failure(request, 'failed', 'Pinned full reference sequence does not match request')
        p = request.identity['residue_1based']
        if not 1 <= p <= len(sequence) or sequence[p-1] != request.identity['ref_aa']:
            return b'', self._failure(request, 'unavailable', 'Exact reference residue correspondence failed')
        mapping = model.get('mapping')
        if mapping == 'fragment':
            try: require_sha(model.get('mapping_qualification_sha256'))
            except ValueError: return b'', self._failure(request, 'unavailable', 'Fragment requires a qualified explicit reference-offset map')
            start = model.get('reference_start_1based')
            if type(start) is not int or start < 1:
                return b'', self._failure(request, 'unavailable', 'Fragment reference offset is missing')
        elif mapping == 'full' and model.get('fragmented_accession') is False: start = 1
        else: return b'', self._failure(request, 'unavailable', 'Full/fragment source mapping is unqualified')
        try:
            with self._lock:
                raw, disposition, execution = self._payload(model, request, store, policy)
            parsed = self._parse(raw, model)
            if parsed['entry_id'] != model['model_id'] or parsed['chain'] != model['chain_id']:
                raise ValueError('Model/chain identity differs from qualified source')
            source_seq = parsed['sequence']; n = len(source_seq)
            if mapping == 'full' and (source_seq != sequence or parsed['db_begin'] != 1 or parsed['db_end'] != len(sequence)):
                raise ValueError('Full sequence/model span correspondence failed')
            if mapping == 'fragment' and sequence[start-1:start-1+n] != source_seq:
                raise ValueError('Explicit fragment/reference sequence correspondence failed')
            i = p - start
            if not 0 <= i < n:
                return raw, self._failure(request, 'unavailable', 'Reference residue is outside the qualified fragment')
            q = np.frombuffer(parsed['plddt'],dtype='<f8')
            xyz = np.frombuffer(parsed['coords'],dtype='<f8').reshape(-1,3)
            values = {'pLDDT':float(q[i]), 'RSA':None, **{name:None for name in GEOMETRY_NAMES}}
            units = {'pLDDT':'pLDDT_0_100', 'RSA':'relative_solvent_accessibility',
                **{name:'count' if 'contacts' in name else 'angstrom' for name in GEOMETRY_NAMES}}
            geometry_reason = 'fragment_is_not_whole_protein_geometry' if mapping == 'fragment' else 'too_short_for_nonlocal_geometry'
            provenance = self._provenance(model)
            if mapping == 'full' and n >= 14:
                try:
                    recipe = self.geometry_recipe()
                    values.update({name:float(value) for name,value in zip(GEOMETRY_NAMES,recipe['descriptors'](xyz)[i])})
                    geometry_reason = 'qualified_reference_geometry'
                    provenance['geometry_sha256'] = GEOMETRY_SHA256
                except (ValueError,FileNotFoundError) as exc:
                    geometry_reason = 'unavailable_qualified_geometry_recipe: ' + str(exc)
            metadata = dict(reference_context_only=True, chain_id=parsed['chain'], model_id=parsed['entry_id'],
                mapping=mapping, reference_start_1based=start, reference_end_1based=start+n-1,
                structure_coverage=n/len(sequence), geometry_status=geometry_reason,
                confidence_interpretation='limited_reference_confidence' if q[i] < 70 else 'reference_model_confidence',
                RSA_status='unavailable_no_qualified_pinned_shrake_rupley_dependency_recipe',
                limits=['Reference geometry does not establish mutation damage, disorder or an interface.',
                    'C-alpha distances are not solvent accessibility.'])
            return raw, EvidenceDecision(request,True,'partial','Qualified reference confidence/context; RSA unavailable',
                values,units,provenance,metadata,cache_disposition=disposition,execution=execution)
        except FileNotFoundError as exc: return b'', self._failure(request,'unavailable',str(exc))
        except PermissionError as exc: return b'', self._failure(request,'pending_authorization',str(exc))
        except BudgetExceeded as exc: return b'', self._failure(request,'unavailable',str(exc))
        except (KeyError,ValueError,TypeError,OSError,UnicodeError,EOFError,IndexError,zlib.error) as exc:
            return b'', self._failure(request,'failed',f'Qualified structure extraction failed: {exc}')

    def obtain(self, request, store, policy):
        if request.provider_id != self.provider_id or request.recipe_sha256 != self.recipe_sha256:
            return self._failure(request,'failed','Request does not bind this provider source/recipe')
        if store is None: return self._derive(request, store, policy)[1]
        try: return store.get_or_compute(request, lambda: self._derive(request, store, policy))
        except (ValueError,OSError) as exc:
            return self._failure(request,'failed',f'Evidence cache integrity failure: {exc}')
