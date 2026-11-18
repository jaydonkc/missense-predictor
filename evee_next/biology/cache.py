"""Closed, selective access to the current two-protein development roster.

The only outcome decoder is admitted_targets, called after release, source and
identity checks. Whole-panel outcome readers and historical target donors are
never imported. Source numeric feature members are projected as byte rows before
NumPy conversion, avoiding numeric decoding of excluded protein rows.
"""
from io import BytesIO
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
from .chemistry import Chemistry, ProteinVariant

STUDY=Path('research/evee-next-2026-11-05')
INPUTS={
 'roles':('effect-supervision-clinical-bridge-v1/role-readiness/output-v1/role-readiness.json','17f83e43e16d78315bf6f5bcc7ddc79879a26d0c3ae8baadc28272c7ee823326'),
 'identity':('effect-readiness/abundance-development-data-v1/output/protein-input.json','30b4ddd382182303c5b1d45f83ac299ce0289c125825a1723ba1d627bace4886'),
 'targets':('effect-supervision-clinical-bridge-v1/teacher-fit-v1/training-inputs.npz','27c9b5a52bab8f415f5a0b1fa15d5d160e4d6396d87801a43d8616b54db8bd80'),
 'saprot':('saprot-effect-context-v1/capture/qualified-features.npz','44d164fdc3a70311ab3930bf278b0b9e37bfba398c58e628ccb39fed21ca304b'),
 'protein':('effect-readiness/effect-protein-learning-v1/execution-v3/operation-v1/bounded-readback/features.npz','f937ae7268788929bbb9d046f7e25a22756a6b4e632c1e58bb76b212d48fefd4'),
 'context':('saprot-effect-context-v1/inputs/baseline.npz','68a52e8399c6e140d975325a115b575fd0dac3dcff07a69f63cf79d5f3b936b9'),
 'stability':('missense-effect-hypothesis-v1/output/panel-b/predictions.npz','db110b96e8e3a7e497b2040f7d0698fe086362ba8aac56d8598233a6b1fca9df'),
}


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def checked(path,sha):
    p=Path(path)
    if not p.is_file() or digest(p)!=sha: raise ValueError('Input hash differs: '+str(p))
    return p


def selected_member(path,key,indices):
    raw_indices=np.asarray(indices)
    if raw_indices.dtype.kind not in 'iu':raise ValueError('Exact integer selection indices required')
    ix=np.asarray(indices,dtype=int)
    if ix.ndim!=1 or (ix<0).any() or len(set(ix.tolist()))!=len(ix):
        raise ValueError('Unique nonnegative selected row indices required')
    with zipfile.ZipFile(path) as z, z.open(key+'.npy') as f:
        version=np.lib.format.read_magic(f)
        read={ (1,0):np.lib.format.read_array_header_1_0,(2,0):np.lib.format.read_array_header_2_0}.get(version)
        if read is None:raise ValueError('Unsupported NPY version')
        shape,fortran,dtype=read(f)
        if fortran or dtype.hasobject or not shape or (len(ix) and ix.max()>=shape[0]):
            raise ValueError('Unsupported selected member')
        width=int(np.prod(shape[1:],dtype=int));rowbytes=width*dtype.itemsize
        parts={};offset=0
        for i in sorted(ix.tolist()):
            skipped=f.read((i-offset)*rowbytes)
            if len(skipped)!=(i-offset)*rowbytes:raise ValueError('Truncated skipped rows')
            raw=f.read(rowbytes)
            if len(raw)!=rowbytes:raise ValueError('Truncated selected row')
            parts[i]=raw;offset=i+1
        return np.frombuffer(b''.join(parts[i] for i in ix),dtype=dtype).copy().reshape((len(ix),)+shape[1:])


def aligned(path,expected_ids,keys):
    with np.load(path,allow_pickle=False) as z:ids=z['case_ids'].astype(str)
    if len(set(ids.tolist()))!=len(ids):raise ValueError('Duplicate source feature IDs')
    lookup={x:i for i,x in enumerate(ids)}
    if not set(expected_ids)<=set(lookup):raise ValueError('Missing admitted source features')
    ix=[lookup[x] for x in expected_ids]
    return {key:selected_member(path,key,ix) for key in keys}


def admitted_targets(archive,ids,genes,indices,released):
    if released is not True or not set(genes)<= {'ASPA','G6PD'}:
        raise ValueError('Current two-protein numeric release required')
    for key,expected in [('case_ids',ids),('genes',genes),('case_indices',indices)]:
        if not np.array_equal(archive[key],expected):raise ValueError('Target identity/order differs: '+key)
    y=archive['targets']
    if y.dtype!=np.float64 or y.shape!=(len(ids),) or not np.isfinite(y).all():
        raise ValueError('Invalid permitted target vector')
    return y


def load_inputs(root: Path,release: dict,*,targets=False):
    paths={k:checked(root/STUDY/name,sha) for k,(name,sha) in INPUTS.items()}
    if release.get('input_hashes')!={k:v[1] for k,v in INPUTS.items()}:
        raise ValueError('Release does not bind current inputs')
    role=json.loads(paths['roles'].read_bytes());identity=json.loads(paths['identity'].read_bytes())
    if identity['labels_present'] is not False or role['status']!='complete_target_free_admission_map':
        raise ValueError('Current role and label-free identity contracts required')
    ix=role['effect_target_projection']['permitted_case_indices']
    if ix!=list(range(15519)):raise ValueError('Permitted row census differs')
    rows=[identity['cases'][i] for i in ix];ids=[r['id'] for r in rows];genes=[r['gene'] for r in rows]
    if len(set(ids))!=len(ids) or {g:genes.count(g) for g in set(genes)}!={'ASPA':5843,'G6PD':9676}:
        raise ValueError('Two-source identity census differs')
    references={}
    for gene in ('ASPA','G6PD'):
        s=role['teacher_sources'][gene]
        if s['new_targets_permitted'] is not True or s['new_permitted_target_case_indices']!=[i for i,r in enumerate(rows) if r['gene']==gene]:
            raise ValueError('Current source permission differs')
        seq=identity['targets'][gene]['sequence']
        if hashlib.sha256(seq.encode()).hexdigest()!=s['reference_sha256']:raise ValueError('Full reference differs')
        references[gene]=dict(sequence=seq,sha256=s['reference_sha256'])
    for r in rows:
        ProteinVariant(r['gene'],r['gene'],r['position0']+1,r['ref_aa'],r['alt_aa'],references[r['gene']]['sha256']).validate(references[r['gene']]['sequence'])
    chem=Chemistry.from_project(root);c=np.array([chem.vector(r['ref_aa'],r['alt_aa']) for r in rows])
    sa=aligned(paths['saprot'],ids,['score_U','score_T','ESMC_LLR'])
    pr=aligned(paths['protein'],ids,['features'])['features']
    ctx=aligned(paths['context'],ids,['pLDDT','RSA','RSA_available'])
    ddg=aligned(paths['stability'],ids,['predicted_ddg'])['predicted_ddg']
    mask=ctx['RSA_available'].astype(bool);rsa=ctx['RSA']
    if not np.isfinite(ctx['pLDDT']).all() or not np.isfinite(rsa[mask]).all() or ((rsa[mask]<0)|(rsa[mask]>1)).any():
        raise ValueError('Invalid qualified reference context')
    context=np.column_stack([ctx['pLDDT']/100,np.where(mask,rsa,0.),mask.astype(float)])
    protein=np.column_stack([pr,sa['ESMC_LLR'],sa['score_U'],sa['score_T']-sa['score_U'],ddg])
    if not np.isfinite(protein).all():raise ValueError('Nonfinite protein features')
    result=dict(rows=rows,ids=np.array(ids),genes=np.array(genes),chemistry=c,context=context,
                protein=protein,references=references,rsa_available=mask,
                native_esmc=sa['ESMC_LLR'],native_saprot=sa['score_T'])
    if targets:
        with np.load(paths['targets'],allow_pickle=False) as z:
            result['targets']=admitted_targets(z,ids,genes,ix,release.get('numeric_targets_authorized'))
    return result
