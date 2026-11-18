"""Explicit hash-bound scoring admission, separate from provider availability."""
import hashlib,json
from pathlib import Path
from .routing import require_sha


def pinned_json(path,pin):
    require_sha(pin);data=Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest()!=pin: raise ValueError('Pinned artifact hash differs')
    def pairs(items):
        d={}
        for k,v in items:
            if k in d:raise ValueError('Duplicate JSON field')
            d[k]=v
        return d
    return json.loads(data,object_pairs_hook=pairs,parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))


class EvidenceRegistry:
    def __init__(self,path,manifest_sha256):
        self.path=Path(path);self.pin=manifest_sha256
        value=pinned_json(self.path,self.pin)
        if set(value)!={'schema','entries'} or value['schema']!='routed-evidence-registry/v1' or not isinstance(value['entries'],list):
            raise ValueError('Qualified evidence registry required')
        self.entries={}
        expected={'artifact_sha256','state','checkpoint','feature_map_id','assembly','qualification_sha256','protein_ids'}
        for entry in value['entries']:
            if set(entry)!=expected: raise ValueError('Registry entry has unexpected fields')
            for key in ('artifact_sha256','feature_map_id','qualification_sha256'): require_sha(entry[key])
            if entry['state'] not in {'diagnostic','offline_candidate','development_supported','confirmed_scoring'}:
                raise ValueError('Unknown evidence admission state')
            if entry['artifact_sha256'] in self.entries:raise ValueError('Duplicate scoring artifact')
            if not isinstance(entry['protein_ids'],list) or any(not isinstance(p,str) for p in entry['protein_ids']):
                raise ValueError('Explicit protein scope required')
            self.entries[entry['artifact_sha256']]=dict(entry)

    def require_confirmed(self,artifact_sha256,checkpoint,feature_map_id,qualification_sha256):
        entry=self.entries.get(artifact_sha256)
        if not entry or entry['state']!='confirmed_scoring':
            raise ValueError('Only explicitly confirmed_scoring artifacts may change routine predictions')
        if (entry['checkpoint'],entry['feature_map_id'],entry['qualification_sha256']) != (checkpoint,feature_map_id,qualification_sha256):
            raise ValueError('Scoring registry checkpoint/feature/qualification differs')
        return dict(entry)
