"""Lock a no-fit ASPA transfer diagnostic before joining outcomes."""
import csv
import gzip
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parent
STRUCTURE = ROOT.parent/'model-3d-understanding-v1/nonlocal-contact-v1'
BUNDLE = ROOT.parent/'named-effect-routing-v1/bundle-v1'
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def lock():
    coords = STRUCTURE/'coordinates/2O4H.cif.gz'
    assert sha(coords)==json.loads((STRUCTURE/'coordinates/2O4H-receipt.json').read_text())['compressed_sha256']
    packet = BUNDLE/'abundance.json'
    manifest = json.loads((BUNDLE/'manifest.json').read_text())
    assert sha(packet)==manifest['endpoints']['abundance']['packet_sha256']
    reference = manifest['references']['ASPA']
    entity = json.loads((STRUCTURE/'sources/entity-2O4H_1.json').read_text())
    sequence = entity['entity_poly']['pdbx_seq_one_letter_code_can']
    assert sequence == reference['sequence']
    aa3 = dict(zip('ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCEQGHILKMFPSTWYV',strict=True))
    # Validate the projected atom table fields before indexing this exact pinned CIF.
    labels = []
    atoms, ligands = defaultdict(list), defaultdict(list)
    with gzip.open(coords,'rt') as stream:
        for line in stream:
            if line.startswith('_atom_site.'):
                labels.append(line.strip().split('.')[1])
            if not line.startswith(('ATOM ','HETATM ')):
                continue
            row = dict(zip(labels,line.split(),strict=True))
            if row['pdbx_PDB_model_num']!='1' or row['auth_asym_id']!='A' or row['type_symbol'] in {'H','D'} or row['label_alt_id'] not in {'.','A'} or float(row['occupancy'])<=0:
                continue
            coord = [float(row['Cartn_'+axis]) for axis in 'xyz']
            if row['group_PDB']=='ATOM' and row['label_entity_id']=='1':
                pos=int(row['label_seq_id'])
                assert sequence[pos-1]==aa3[row['label_comp_id']]
                atoms[pos].append(coord)
            if row['label_comp_id'] in {'ZN','AS9'}:
                ligands[row['label_comp_id']].append(coord)
    assert len(ligands['ZN'])==1 and len(ligands['AS9'])==13
    all_ligands=ligands['ZN']+ligands['AS9']
    distances={pos: float(cdist(points,all_ligands).min()) for pos,points in atoms.items()}
    data=io.StringIO()
    writer=csv.writer(data,delimiter='\t',lineterminator='\n')
    writer.writerow(['substitution','abundance_percentile','contact_distance_A','category','abundance_only'])
    for gene,pos,ref,alt,score,_ in json.loads(packet.read_text()):
        if gene!='ASPA':
            continue
        assert sequence[pos-1]==ref
        dist=distances.get(pos)
        base='dual' if score<0.5 else 'neither'
        category=base if dist is not None else 'unavailable'
        if dist is not None and score>=0.5 and dist<=5:
            category='selective'
        writer.writerow((f'{ref}{pos}{alt}',score,'' if dist is None else dist,category,base))
    raw=data.getvalue().encode()
    payload=gzip.compress(raw,mtime=0)
    path=ROOT/'predictions.tsv.gz'
    with path.open('xb') as fp: fp.write(payload)
    receipt=dict(source_structure_sha256=sha(coords),structured_residues=len(atoms),target_free_packet=sha(packet),
        predictions_sha256=sha(path),reference=reference['protein_id'],protein_sha256=reference['sha256'],
        script_sha256=sha(Path(__file__)),protocol_sha256=sha(ROOT/'PROTOCOL.md'),rows=len(raw.splitlines())-1)
    with (ROOT/'score-lock.json').open('x') as fp: json.dump(receipt,fp,indent=2);fp.write('\n')
if __name__=='__main__':
    lock()
