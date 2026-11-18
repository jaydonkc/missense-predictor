"""Target-free receiving math for one explicitly chosen retained checkpoint.

No backbone, cohort loader, fitting, reference retrieval or outcome access occurs
here. The manifest SHA is supplied by the caller. A model bundle binds model
state; the caller separately binds every query's genomic/protein/source identity
and qualifies its raw extraction before using these coordinates in a report.
"""
import ast
from io import BytesIO
import hashlib
import json
import math
from pathlib import Path
import struct
from types import MappingProxyType
import zipfile

import numpy as np

from .routing import BaselinePrediction, BaselineResult, BaselineState, require_sha


CUTOFF = .2787385170622812
_STUDY = 'research/evee-next-2026-11-05/'
_SOURCES = {
    'reader': (_STUDY+'dna-first-moment-v1/checkpoint_f32.py',
        'b38ecf0d842f515b14ada2bf1b0015ff02ba1ceead4b4bc1f3949b308c20ce59'),
    'means': (_STUDY+'dna-first-moment-v1/mean_features.py',
        'd27a24e0c713640ccbee52e645a3abad796d87fd6460654226b5e956a8dbb5d9'),
    'native': (_STUDY+'strong-dna-fusion-train-v1/dna_core.py',
        'd3952a510c555979f6e1e5516b281a5342879cb79a4c06a8e5321107316cb751'),
    'protein': (_STUDY+'strong-dna-fusion-train-v1/model.py',
        '99636b85cb9f20c34e0376bdd4408b289ec370db84b98d48f7a7fa487002bedc'),
    'chemistry': ('evee_next/cached.py',
        'a7db8ab02bc081bdc55e53ac01bef22fae0b31f3aa1ebfbe885a7cf59a5f8b73'),
}
_SHAPES = dict(protein_mean=(512,), protein_scale=(512,),
    base_median=(5324,), base_mean=(5324,), base_scale=(5324,),
    scalar_mean=(9,), scalar_scale=(9,), W0=(5333,64), b0=(64,))
_STATE_NAMES = set(_SHAPES) | {'W1','b1'}
_FILES = {'state.npz','dna_projection.safetensors','chemistry_tables.json'}


def _require(ok, why='Retained serving-state contract differs'):
    if not ok:
        raise ValueError(why)


def _frozen(value):
    value = np.asarray(value)
    return np.frombuffer(value.tobytes(order='C'), dtype=value.dtype).reshape(value.shape)


def _vector(value, width, name, *, missing=False):
    a = np.asarray(value)
    _require(a.dtype.kind == 'f' or a.dtype.kind in 'iu', name+' must be real numeric coordinates')
    a = np.asarray(a, dtype=np.float64)
    _require(a.shape == (width,) and not np.isinf(a).any()
             and (missing or np.isfinite(a).all()), name+' shape/finite values differ')
    return a


def _source_functions(key, names, namespace=None):
    """Compile named definitions only; none of a study module's imports execute."""
    relative, pin = _SOURCES[key]
    path = Path(__file__).resolve().parents[2] / relative
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ValueError('Pinned serving math source unavailable: '+relative) from exc
    _require(hashlib.sha256(raw).hexdigest() == pin, 'Pinned serving math source changed: '+relative)
    tree = ast.parse(raw, filename=str(path))
    definitions = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name in names]
    _require({n.name for n in definitions} == set(names), 'Pinned definitions incomplete')
    environment = dict(np=np, Path=Path, json=json, math=math, struct=struct, require=_require)
    environment.update(namespace or {})
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(path), 'exec'), environment)
    return environment


def _json(raw):
    def unique(pairs):
        result = {}
        for k,v in pairs:
            _require(k not in result, 'Duplicate JSON member')
            result[k] = v
        return result
    try:
        return json.loads(raw, object_pairs_hook=unique,
                          parse_constant=lambda _: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError('Invalid serving-state JSON') from exc


class ServingState:
    """One original fold's frozen model and scalers, trusted by a manifest SHA.

    ``transform`` returns a read-only (5333,) float64 vector. ``score`` returns a
    float sigmoid score. ``prediction`` additionally creates private routed state;
    it does not establish query identity or extraction qualification.
    """
    def __init__(self, bundle: Path, manifest_sha256: str):
        require_sha(manifest_sha256)
        self.bundle = Path(bundle)
        manifest_path = self.bundle/'manifest.json'
        _require(not manifest_path.is_symlink() and manifest_path.stat().st_size <= 2**20,
                 'Bounded regular serving manifest required')
        raw = manifest_path.read_bytes()
        _require(hashlib.sha256(raw).hexdigest() == manifest_sha256, 'Serving manifest hash differs')
        m = _json(raw)
        _require(isinstance(m,dict) and set(m) == {'schema','fold','feature_map_id','checkpoint_sha256',
                 'files','source_bindings','recipes'} and m['schema'] == 'retained-serving-state/v1')
        _require(type(m['fold']) is int and m['fold'] in (0,1,2), 'One original fold must be explicit')
        for name in ('feature_map_id','checkpoint_sha256'):
            require_sha(m[name])
        for name in ('source_bindings','recipes'):
            _require(isinstance(m[name],dict) and bool(m[name]), 'Nonempty '+name+' hashes required')
            for k,v in m[name].items():
                _require(isinstance(k,str) and bool(k), 'Named source/recipe required')
                require_sha(v)
        _require(isinstance(m['files'],dict) and set(m['files']) == _FILES, 'Only target-free serving files allowed')
        payload = {}
        for name, cap in [('state.npz',16*2**20),('dna_projection.safetensors',8*2**20),
                          ('chemistry_tables.json',2**20)]:
            spec = m['files'][name]
            _require(isinstance(spec,dict) and set(spec) == {'sha256','bytes'}, 'Exact file byte/hash specification required')
            require_sha(spec['sha256'])
            _require(type(spec['bytes']) is int and 0 < spec['bytes'] <= cap, 'Serving file size exceeds bound')
            path = self.bundle/name
            _require(not path.is_symlink() and path.is_file() and path.stat().st_size == spec['bytes'],
                     'Serving file identity/byte count differs: '+name)
            payload[name] = path.read_bytes()
            _require(hashlib.sha256(payload[name]).hexdigest() == spec['sha256'], 'Serving file hash differs: '+name)
        with zipfile.ZipFile(BytesIO(payload['state.npz'])) as z:
            members = z.infolist()
            _require(len(members) == len(_STATE_NAMES) and {n.filename for n in members} == {k+'.npy' for k in _STATE_NAMES}
                     and sum(n.file_size for n in members) <= 8*2**20, 'Only bounded named state arrays are allowed')
        with np.load(BytesIO(payload['state.npz']),allow_pickle=False) as z:
            arrays = {k:z[k] for k in _STATE_NAMES}
        for k,a in arrays.items():
            shape_ok = a.shape == _SHAPES[k] if k in _SHAPES else (
                a.shape in ((64,), (64,1)) if k == 'W1' else a.shape in ((),(1,)))
            _require(a.dtype == np.dtype('float64') and shape_ok and np.isfinite(a).all(),
                     'Finite float64 serving-state shape required: '+k)
            if k.endswith('scale'):
                _require(np.all(a > 0), 'Saved positive scales required: '+k)
        arrays['W1'] = arrays['W1'].reshape(64)
        arrays['b1'] = arrays['b1'].reshape(())
        h = hashlib.sha256()
        for k in ('W0','b0','W1','b1'):
            h.update(str(arrays[k].shape).encode()); h.update(np.ascontiguousarray(arrays[k]).tobytes())
        _require(h.hexdigest() == m['checkpoint_sha256'], 'Actual neural parameter identity differs')
        self.state = MappingProxyType({k:_frozen(a) for k,a in arrays.items()})
        reader = _source_functions('reader', {'unique_object','load_file'})
        # Recheck bytes after the pure reader's bounded path access.
        projection_path = self.bundle/'dna_projection.safetensors'
        p = reader['load_file'](projection_path)
        _require(projection_path.read_bytes() == payload['dna_projection.safetensors'], 'Projection changed during read')
        expected = {f'{side}.{d}.{kind}': ((64,4096) if kind == 'weight' else (64,))
                    for d in range(2) for side in ('left','right') for kind in ('weight','bias')}
        expected.update(epsilon_eye=(64,64), **{'readout.weight':(2,4608),'readout.bias':(2,)})
        _require(set(p) == set(expected) and all(p[k].shape == shape for k,shape in expected.items()),
                 'Complete unconditioned late projection state required')
        self.raw_projection_state = MappingProxyType({k:_frozen(a) for k,a in p.items()})
        self.tables = _json(payload['chemistry_tables.json'])
        self._means = _source_functions('means', {'dna_array','mean_delta','weights','project_means'})
        self._chemistry = _source_functions('chemistry', {'chemistry_vector'}, {'AA':'ACDEFGHIKLMNPQRSTVWY'})['chemistry_vector']
        self.fold = m['fold']
        self.feature_map_id = m['feature_map_id']
        self.checkpoint_sha256 = m['checkpoint_sha256']
        self.provenance = MappingProxyType(dict(manifest_sha256=manifest_sha256, fold=self.fold,
            feature_map_id=self.feature_map_id, checkpoint_sha256=self.checkpoint_sha256,
            source_bindings=MappingProxyType(dict(m['source_bindings'])), recipes=MappingProxyType(dict(m['recipes'])),
            math_sources=MappingProxyType({path:pin for path,pin in _SOURCES.values()}),
            scope='explicit_checkpoint_development', cutoff=CUTOFF))

    def chemistry(self, ref, alt):
        """Original 456-coordinate canonical substitution descriptor."""
        _require(isinstance(ref,str) and isinstance(alt,str)
                 and ref in 'ACDEFGHIKLMNPQRSTVWY' and alt in 'ACDEFGHIKLMNPQRSTVWY'
                 and len(ref) == len(alt) == 1 and ref != alt, 'Canonical nonidentical residues required')
        try:
            value = self._chemistry(ref,alt,self.tables)
        except (KeyError,TypeError,ValueError) as exc:
            raise ValueError('Frozen chemistry table is incomplete') from exc
        return _frozen(_vector(value,456,'Chemistry'))

    def signed_means(self, raw_dna):
        """Float64 selected ALT-minus-REF mean projections; no projection bias."""
        x = self._raw_dna(raw_dna)
        means = self._means['mean_delta'](x)
        return _frozen(self._means['project_means'](means,self.raw_projection_state)[0])

    @staticmethod
    def _raw_dna(raw_dna):
        x = np.asarray(raw_dna)
        _require(x.shape == (1,2,2,256,4096) and x.dtype.kind == 'f' and np.isfinite(x).all(),
                 'One decoded finite floating-point 2x2x256x4096 DNA capture required')
        return x

    def native_representations(self, raw_dna, raw_protein512, *, device='cpu'):
        """Pinned original Torch float32 DNA/protein math; no model inference.

        Returns (D4096, P512, M256). Torch is imported only here. Missing Torch
        fails explicitly; no installation or alternate precision is selected.
        """
        x = self._raw_dna(raw_dna)
        p = _vector(raw_protein512,512,'Raw projected protein')
        try:
            import torch
        except ImportError as exc:
            raise RuntimeError('Torch is required for original native float32 representation math') from exc
        _require(torch.backends.cuda.matmul.allow_tf32 is False
                 and torch.backends.cudnn.allow_tf32 is False
                 and torch.get_float32_matmul_precision() == 'highest',
                 'Original native projection precision requires TF32 disabled and highest float32 precision')
        env = _source_functions('native', {'finite','matrix_transform','ReferenceFusion'},
                                {'torch':torch,'nn':torch.nn})
        env = _source_functions('protein', {'protein_view','VariantModel'},
            {'torch':torch,'nn':torch.nn,'F':torch.nn.functional,
             'finite':env['finite'],'ReferenceFusion':env['ReferenceFusion']})
        _require(str(device) == 'cpu' or str(device).startswith('cuda'), 'Explicit CPU/CUDA device required')
        # Initialization is immediately replaced by saved state; preserve RNG.
        devices = [] if str(device) == 'cpu' else [torch.device(device).index or torch.cuda.current_device()]
        with torch.random.fork_rng(devices=devices):
            model = env['VariantModel']('late').to(device)
        model.load_state_dict({k:torch.tensor(np.array(v,copy=True),device=device)
                               for k,v in self.raw_projection_state.items()},strict=True)
        model.eval().requires_grad_(False)
        q = ((p-self.state['protein_mean'])/self.state['protein_scale']).astype(np.float32)
        with torch.inference_mode():
            result = model.features(torch.as_tensor(np.asarray(x,dtype=np.float32),device=device),
                                    torch.as_tensor(q[None],device=device)).cpu().numpy()[0]
        _require(result.shape == (4608,) and result.dtype == np.float32 and np.isfinite(result).all(),
                 'Original native representation shape/precision differs')
        return _frozen(result[:4096]), _frozen(result[4096:]), self.signed_means(x)

    def transform(self, raw_dna, raw_protein512, chemistry456, Q2, S2, masked_llr, *, device='cpu'):
        """Apply original native math and all frozen transforms to one query."""
        D,P,M = self.native_representations(raw_dna,raw_protein512,device=device)
        return self.transform_from_representations(D,P,M,chemistry456,Q2,S2,masked_llr)

    def transform_from_representations(self, D4096, P512, M256, C456, Q2, S2, LLR):
        """Assemble original representations in fixed order and apply saved scales.

        P512 must already include original protein standardization plus native
        float32 LayerNorm. D4096 must be the native late state representation.
        This receiver does not certify either upstream extraction.
        """
        d,p,m,c,q,s = [_vector(v,w,n,missing=n=='Structure') for v,w,n in [
            (D4096,4096,'DNA'),(P512,512,'Protein'),(M256,256,'Means'),
            (C456,456,'Chemistry'),(Q2,2,'Availability'),(S2,2,'Structure')]]
        _require(set(q) <= {0.,1.} and q[1] <= q[0], 'Binary structure/burial availability required')
        known = q[1] == 1
        _require((np.isfinite(s).all() if known else np.isnan(s).all()),
                 'Burial availability and structural missingness differ')
        if known:
            _require(s[0] in (0.,1.) and s[1] >= 0 and (s[0] == 1) == (s[1] < .2),
                     'Known reference burial/RSA differs')
        _require(np.array_equal(c[443:450:2],np.abs(c[442:450:2])), 'Signed/absolute chemistry delta order differs')
        llr = np.asarray(LLR)
        _require(llr.dtype.kind in 'fiu' and llr.shape == ()
                 and np.isfinite(llr), 'Finite scalar masked LLR required')
        interaction = (1-np.clip(s[1],0,1))*c[442:450] if known else np.zeros(8)
        base = np.concatenate((d,p,c,q,s,m))
        base = (np.where(np.isnan(base),self.state['base_median'],base)-self.state['base_mean'])/self.state['base_scale']
        scalars = (np.r_[float(LLR),interaction]-self.state['scalar_mean'])/self.state['scalar_scale']
        result = np.concatenate((base,scalars))
        _require(result.shape == (5333,) and np.isfinite(result).all(), 'Complete finite transformed coordinates required')
        return _frozen(result)

    def score(self, X):
        """Original float64 5333→64 ReLU→1 sigmoid replay, returning one float."""
        x = _vector(X,5333,'Transformed input')
        hidden = np.maximum(0,x@self.state['W0']+self.state['b0'])
        logit = hidden@self.state['W1']+self.state['b1']
        _require(np.isfinite(hidden).all() and np.isfinite(logit), 'Nonfinite frozen neural replay')
        return float(np.exp(-np.logaddexp(0,-logit)))

    def prediction(self, X, variant_key):
        """Create routed private state after the caller validates query identity."""
        _require(isinstance(variant_key,str) and bool(variant_key), 'Bound query identity required')
        x = _vector(X,5333,'Transformed input')
        prediction = BaselinePrediction('scored',self.score(x),CUTOFF,self.checkpoint_sha256,
            self.feature_map_id,'explicit_checkpoint_development',
            reason='retained original checkpoint; fresh source/extraction qualification is supplied separately')
        private = BaselineState(variant_key,self.checkpoint_sha256,self.feature_map_id,
            tuple(x@self.state['W0']+self.state['b0']),tuple(self.state['W1']),float(self.state['b1']))
        return BaselineResult(prediction,private)
