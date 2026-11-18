"""Synthetic data for the runnable preview."""

from dataclasses import asdict
import hashlib
import json

import numpy as np

from evee_next.biology.chemistry import Chemistry
from evee_next.report import MissenseMapping, VariantRequest, normalize

MAPPING = MissenseMapping('ENST00001.1', 'ENSP00001.1', '+', 2, 'A', 'D')

def chemistry():
    tables = {name: {"A": 1., "D": 2., "E": 3.}
              for name in ("charge", "volume", "polarity", "hydrophobicity")}
    tables.update(grantham={"A": {"D": 126.}, "E": {"A": 107.}},
                  pam250={"A": {"D": 0.}, "E": {"A": 0.}})
    return Chemistry(tables, "synthetic_chemistry")

def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()

def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n')
    return hashlib.sha256(path.read_bytes()).hexdigest()

def fixture(directory):
    """Literal retained-width vectors whose first two hidden units give logit3.3."""
    seq = 'MAAAA'
    refsha = hashlib.sha256(seq.encode()).hexdigest()
    rows = []
    folds = []
    for fold, (pos, ref, alt) in enumerate([(10,'A','G'), (11,'C','T'), (12,'G','A')]):
        variant = normalize(VariantRequest('1', pos, ref, alt))
        rows.append(dict(source_row=100+fold, panel_row=fold, fold=fold,
                         group='group'+str(fold), component='component'+str(fold),
                         variant=asdict(variant), mapping=asdict(MAPPING),
                         reference_sha256=refsha))
        X = np.zeros((1,5333)); X[0,[0,4096,5068]] = [2.,-1.,1.]
        W0 = np.zeros((5333,64)); W0[0,:2] = [1.,-1.]
        W0[4096,:2] = [2.,2.]; W0[5068,:2] = [1.,0.]
        b0 = np.zeros(64); b0[:2] = [.5,-.5]
        W1 = np.zeros((64,1)); W1[:2,0] = [2.,-3.]
        archive = directory/f'fold{fold}.npz'
        np.savez_compressed(archive, X=X, W0=W0, b0=b0, W1=W1, b1=np.array([.3]))
        feature_map = dict(fold=fold, original_fitting_matrix_sha256=str(fold+1)*64,
                           original_heldout_matrix_sha256=str(fold+4)*64,
                           input_dimension=5333, original_neural_source_pin='a'*64,
                           coordinate_semantics='original_fitting_only_transformed_coordinates')
        folds.append(dict(fold=fold, feature_map=feature_map,
                          feature_map_id=hashlib.sha256(canonical(feature_map)).hexdigest(),
                          artifact=dict(path=archive.name, sha256=hashlib.sha256(archive.read_bytes()).hexdigest(), bytes=archive.stat().st_size)))
    manifest = dict(schema='biology-retained-neural-cache/v1',
                    evaluation_scope='cached_original_heldout_TRAIN_development',
                    target_free=True, cutoff=.5, rows=rows, folds=folds,
                    references={'ENSP00001.1':dict(sequence=seq, sha256=refsha)},
                    chemistry_tables=chemistry().tables,
                    source_provenance={'identity':'synthetic_archived_qualified_mapping_no_fresh_remap'})
    return manifest, write_json(directory/'manifest.json', manifest)
