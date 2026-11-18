"""Target-free access to a fixed roster in the retained neural coordinate maps.

This reader does no feature extraction, transformation fitting, or remapping.
Each cached allele stays with its original held-out TRAIN head. An explicit
manifest SHA is the trust anchor. Original input/transform fingerprints and
reference sequences are part of that manifest; outcome files are not dependencies.
"""
import argparse
from dataclasses import asdict
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
from types import MappingProxyType
import zipfile

import numpy as np

from evee_next.interpretation.neural import FeatureBlock, FrozenReluHead, immutable
from evee_next.report import IdentityResolution, MissenseMapping, NormalizedVariant, VariantRequest, normalize
from .chemistry import Chemistry
from .neural_pipeline import NeuralBiologyPipeline, NeuralFeatureRecord


BLOCKS = (
    FeatureBlock('DNA_bilinear',0,4096,'task_trained_DNA'),
    FeatureBlock('protein_embedding',4096,4608,'biologically_unresolved_embedding'),
    FeatureBlock('amino_acid_chemistry',4608,5064,'amino_acid_chemistry'),
    FeatureBlock('availability',5064,5066,'availability_indicator'),
    FeatureBlock('reference_structure',5066,5068,'reference_structure_context'),
    FeatureBlock('DNA_means',5068,5324,'task_trained_DNA'),
    FeatureBlock('sequence_likelihood',5324,5325,'sequence_likelihood_not_mechanism'),
    FeatureBlock('reference_context_x_chemistry',5325,5333,'reference_context_x_chemistry'),
)
SOURCES = MappingProxyType(dict(zip((b.name for b in BLOCKS),(
    'DNA','direct_protein','amino_acid_chemistry','availability',
    'reference_structure','DNA','sequence_likelihood','reference_context_x_chemistry'))))
SCHEMA = 'biology-retained-neural-cache/v1'
SCOPE = 'cached_original_heldout_TRAIN_development'


def require(ok, message):
    if not ok:raise ValueError(message)


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',', ':'),allow_nan=False).encode()).hexdigest()


def _is_sha(value):
    return isinstance(value,str) and re.fullmatch('[0-9a-f]{64}',value) is not None


def _payload(path, expected):
    require(_is_sha(expected),'A manifest/input SHA256 trust anchor is required')
    try:payload=Path(path).read_bytes()
    except OSError as error:raise ValueError('Cached retained input unavailable: '+str(path)) from error
    require(hashlib.sha256(payload).hexdigest()==expected,'Cached retained input hash differs: '+str(path))
    return payload


def _manifest(bundle, pin):
    value=json.loads(_payload(Path(bundle)/'manifest.json',pin))
    require(set(value)=={'schema','evaluation_scope','target_free','cutoff','rows','folds','references','chemistry_tables','source_provenance'},'Unexpected retained manifest members')
    require(value['schema']==SCHEMA and value['evaluation_scope']==SCOPE
            and value['target_free'] is True and type(value['cutoff']) in (int,float)
            and np.isfinite(value['cutoff']) and 0<=value['cutoff']<=1,'Target-free cached development manifest required')
    return value


class CachedRetainedProvider:
    """Exact admitted cached identities and already transformed 5,333 inputs.

    One instance binds one head/fold. The class also provides the cached identity
    interface for NeuralBiologyPipeline. The archived mapping is qualified by its
    original sources; it is not a fresh genome/codon remap or novel inference.
    Missing cached alleles remain unscored; no other head supplies a substitute.
    """
    def __init__(self,bundle: Path,manifest_sha256: str,fold: int):
        require(type(fold) is int and fold in (0,1,2),'An original retained fold is required')
        self.bundle=Path(bundle)
        manifest=_manifest(self.bundle,manifest_sha256)
        require(isinstance(manifest['rows'],list) and bool(manifest['rows'])
                and isinstance(manifest['folds'],list) and len(manifest['folds'])==3,'Retained identities and three original heads required')
        references=manifest['references'];require(isinstance(references,dict) and references,'Exact reference sequences required')
        for protein,ref in references.items():
            require(re.fullmatch(r'ENSP\d+\.\d+',protein) is not None and set(ref)=={'sequence','sha256'}
                    and isinstance(ref['sequence'],str) and bool(ref['sequence'])
                    and hashlib.sha256(ref['sequence'].encode()).hexdigest()==ref['sha256'],'Reference sequence hash/version differs')
        by_fold={}; map_ids=set()
        for record in manifest['folds']:
            require(set(record)=={'fold','feature_map','feature_map_id','artifact'},'Unexpected retained head fields')
            f=record['fold']; require(type(f) is int and f in (0,1,2) and f not in by_fold,'Duplicate/invalid retained head fold')
            fm=record['feature_map'];require(isinstance(fm,dict) and fm.get('fold')==f
                and fm.get('input_dimension')==5333
                and fm.get('coordinate_semantics')=='original_fitting_only_transformed_coordinates'
                and all(_is_sha(fm.get(k)) for k in ('original_fitting_matrix_sha256','original_heldout_matrix_sha256','original_neural_source_pin'))
                and canonical_sha(fm)==record['feature_map_id'] and record['feature_map_id'] not in map_ids,'Original fold/feature-map binding differs')
            artifact=record['artifact']; require(set(artifact)=={'path','sha256','bytes'}
                and artifact['path']==f'fold{f}.npz' and _is_sha(artifact['sha256'])
                and type(artifact['bytes']) is int and 0<artifact['bytes']<64*1024**2,'Retained artifact path/hash binding differs')
            by_fold[f]=record;map_ids.add(record['feature_map_id'])
        require(set(by_fold)=={0,1,2},'Missing original retained head')
        rows=[]; seen={k:set() for k in ('variant_mapping','source','panel')}; group_folds={};component_folds={};allele_folds={}
        for row in manifest['rows']:
            require(set(row)=={'source_row','panel_row','fold','group','component','variant','mapping','reference_sha256'},'Unexpected cached identity fields')
            f=row['fold'];require(type(f) is int and f in by_fold,'Cached held-out fold differs')
            for k in ('source_row','panel_row'):require(type(row[k]) is int and row[k]>=0,'Invalid original/panel row')
            require(isinstance(row['group'],str) and row['group'] and isinstance(row['component'],str) and row['component'],'Original group/current component required')
            variant=NormalizedVariant(**row['variant']);mapping=MissenseMapping(**row['mapping'])
            require(normalize(VariantRequest(variant.chromosome,variant.position_1based,variant.ref,variant.alt,assembly=variant.assembly))==variant,'Normalized allele metadata differs')
            ref=references.get(mapping.protein_id)
            require(ref is not None and ref['sha256']==row['reference_sha256']
                    and 0<mapping.residue_1based<=len(ref['sequence'])
                    and ref['sequence'][mapping.residue_1based-1]==mapping.ref_aa,'Cached protein/reference/focal residue differs')
            for key,item in [('variant_mapping',(variant,mapping)),('source',row['source_row']),('panel',row['panel_row'])]:
                require(item not in seen[key],'Duplicate cached allele/source/panel identity');seen[key].add(item)
            for lookup,key in [(group_folds,row['group']),(component_folds,row['component']),(allele_folds,variant.variant_id)]:
                require(lookup.setdefault(key,f)==f,'Group/component spans original held-out folds')
            if f==fold:rows.append((row,variant,mapping))
        require(rows and [r[0]['panel_row'] for r in rows]==sorted(r[0]['panel_row'] for r in rows),'Nonempty ordered original held-out rows required')
        record=by_fold[fold];artifact=record['artifact'];payload=_payload(self.bundle/artifact['path'],artifact['sha256'])
        require(len(payload)==artifact['bytes'],'Retained artifact byte count differs')
        allowed={'X','W0','b0','W1','b1'}
        with zipfile.ZipFile(BytesIO(payload)) as archive:
            names=archive.namelist();require(len(names)==5 and set(names)=={k+'.npy' for k in allowed},'Unexpected/duplicate retained archive member or outcomes')
        with np.load(BytesIO(payload),allow_pickle=False) as archive:
            arrays={k:immutable(archive[k]) for k in allowed}
        n=len(rows)
        require(arrays['X'].shape==(n,5333) and arrays['W0'].shape==(5333,64)
                and arrays['b0'].shape==(64,) and arrays['W1'].shape==(64,1)
                and arrays['b1'].shape==(1,),'Retained cached input/head dimensions differ')
        self.head=FrozenReluHead(arrays['W0'],arrays['b0'],arrays['W1'][:,0],float(arrays['b1'][0]),manifest['cutoff'],BLOCKS,record['feature_map_id'])
        self._features=arrays['X'];self._lookup={}
        for i,(row,variant,mapping) in enumerate(rows):
            self._lookup.setdefault(variant.variant_id,[]).append((i,row,variant,mapping))
        self.fold=fold;self.manifest_sha256=manifest_sha256
        self.references=MappingProxyType({k:MappingProxyType(dict(v)) for k,v in references.items()})
        self.chemistry=Chemistry(manifest['chemistry_tables'],'retained_target_free_manifest#sha256='+manifest_sha256)
        self.provenance=dict(manifest['source_provenance'],retained_manifest_sha256=manifest_sha256,
                             identity_scope='archived_qualified_cached_mapping_no_fresh_genomic_or_codon_remap',
                             inference_scope=SCOPE,original_heldout_fold=str(fold))

    def validate_reference(self,protein_id: str,sequence: str):
        ref=self.references.get(protein_id)
        require(ref is not None and sequence==ref['sequence']
                and hashlib.sha256(sequence.encode()).hexdigest()==ref['sha256'],'Protein/reference sequence differs')

    def resolve_identity(self,variant: NormalizedVariant):
        items=self._lookup.get(variant.variant_id)
        if items is None:
            return IdentityResolution((),'outside_cached_original_heldout_roster','outside_admitted_roster',('cached_retained_inputs',))
        require(all(variant==item[2] for item in items),'Cached normalized genomic identity differs')
        return IdentityResolution(tuple(item[3] for item in items),'archived_qualified_cached_missense_no_fresh_remap',
            unresolved_questions=('Only archived cached transcript mappings are represented; other transcript consequences were not assessed here.',))

    def mappings(self,variant: NormalizedVariant):
        return self.resolve_identity(variant).mappings

    def features(self,variant: NormalizedVariant,mapping: MissenseMapping,required_blocks: tuple[str,...]):
        require(required_blocks==tuple(SOURCES),'Required retained block order differs')
        items=self._lookup.get(variant.variant_id)
        if items is None:
            return NeuralFeatureRecord(variant,mapping,self.head.feature_map_id,{},SOURCES,self.provenance)
        matches=[item for item in items if item[2]==variant and item[3]==mapping]
        require(len(matches)==1,'Cached genomic/protein mapping differs')
        i,row,original,original_mapping=matches[0]
        require(variant==original and mapping==original_mapping,'Cached genomic/protein mapping differs')
        provenance=dict(self.provenance,original_source_row=str(row['source_row']),original_panel_row=str(row['panel_row']),
                        reference_sha256=row['reference_sha256'],original_group=row['group'],current_component=row['component'])
        blocks={b.name:self._features[i,b.start:b.end] for b in BLOCKS}
        return NeuralFeatureRecord(variant,mapping,self.head.feature_map_id,blocks,SOURCES,provenance)


def main():
    p=argparse.ArgumentParser(description='Target-free retained-head report for its admitted cached TRAIN roster')
    p.add_argument('bundle',type=Path);p.add_argument('manifest_sha256')
    p.add_argument('chromosome');p.add_argument('position_1based',type=int);p.add_argument('ref');p.add_argument('alt')
    p.add_argument('--transcript-id');p.add_argument('--protein-id')
    args=p.parse_args();request=VariantRequest(args.chromosome,args.position_1based,args.ref,args.alt,args.transcript_id,args.protein_id)
    manifest=_manifest(args.bundle,args.manifest_sha256)
    try:variant=normalize(request)
    except (ValueError,NotImplementedError):variant=None
    matches=[r for r in manifest['rows'] if variant is not None and r['variant']['variant_id']==variant.variant_id]
    require(len({r['fold'] for r in matches})<=1,'Cached allele mappings span incompatible original heads')
    fold=matches[0]['fold'] if matches else 0
    provider=CachedRetainedProvider(args.bundle,args.manifest_sha256,fold)
    report=NeuralBiologyPipeline(provider,provider.chemistry,head=provider.head,feature_provider=provider).analyze(request)
    print(json.dumps(asdict(report),sort_keys=True,indent=2,allow_nan=False))


if __name__=='__main__':main()
