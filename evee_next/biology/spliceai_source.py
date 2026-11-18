"""Pinned, outcome-free SpliceAI VCF/TBI/BGZF readers.

The range mechanics derive from the historical retrieve_scores.py only. This
module has no cohort discovery, DEV filter, target join or promotion rule.
Byte ranges passed to BoundedFetcher are inclusive, as required by HTTP.
"""
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import struct
from types import MappingProxyType
import zlib

from .routing import require_sha, _frozen

SOURCE_SCHEMA = 'spliceai-source/v1'
MANIFEST_SCHEMA = {
    'schema': SOURCE_SCHEMA, 'required': ('source_id', 'assembly', 'reference_sha256',
    'annotation_sha256', 'reference_hash_key', 'release', 'model', 'annotation', 'masked', 'distance',
    'coordinate_system', 'source_scope', 'mane_records'),
    'local': ('score_sha256',),
    'indexed': ('index_url', 'index_sha256', 'score_url', 'score_etag', 'score_bytes'),
}
# Historical source pins are documented identifiers, never implicit approval or
# a configured source. A caller must supply a SHA-pinned complete manifest.
ENSEMBL_MANE14_PINS = MappingProxyType({
    'score_url': 'https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz',
    'score_etag': '"6ab41f97c-64146427056f6"',
    'score_bytes': 28643031420,
    'index_sha256': 'a64c62434c554f5c98a10f84dab8e5544f45b145cd28acf10411ec350e49c20d',
    'mane_gzip_sha256': '4b3992457556e302a5e47e18a305bb4763718377696ca37d5a6b34df7db630d2',
    'header_bgzf_sha256': '7d01fa77eaef0f673c2c847a2bf58e488ce1d6bc06480c0b4b0885d3e3f5e799',
})


def digest(data): return hashlib.sha256(data).hexdigest()

def _require(ok, reason):
    if not ok: raise ValueError(reason)

def _gunzip(data, limit):
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
        result = stream.read(limit + 1)
    _require(len(result) <= limit, 'Decompressed source exceeds bound')
    return result


@dataclass(frozen=True)
class SpliceSourceManifest:
    data: dict
    sha256: str

    def __post_init__(self):
        require_sha(self.sha256)
        d = dict(self.data)
        _require(all(k in d for k in MANIFEST_SCHEMA['required']), 'Incomplete source manifest')
        _require(d['schema'] == SOURCE_SCHEMA, 'Unsupported source schema')
        _require(d['assembly'] == 'GRCh38' and d['coordinate_system'] == 'VCF-1-based', 'Unsupported coordinate/build contract')
        _require(d['reference_hash_key'] in ('genome','reference_sha256','chromosome_fasta_sha256','genome_manifest_sha256'), 'Explicit supported reference hash kind required')
        _require(d['source_scope'] == 'gene', 'This parser only qualifies gene-scope source scores')
        _require(d['annotation'] == 'MANE Select 1.4', 'Unqualified annotation release')
        for k in ('reference_sha256', 'annotation_sha256'): require_sha(d[k])
        for k in ('source_id','release','model'):
            _require(isinstance(d[k],str) and bool(d[k]), 'Missing source/model version')
        _require(type(d['masked']) is bool and type(d['distance']) is int and d['distance'] > 0, 'Explicit mask/distance required')
        for k in ('score_sha256','index_sha256','mane_gzip_sha256','header_bgzf_sha256'):
            if k in d: require_sha(d[k])
        _require(isinstance(d['mane_records'],(list,tuple)), 'MANE records required')
        seen = set()
        for r in d['mane_records']:
            _require(set(r) >= {'transcript_id','gene_id','gene_symbol','strand'}, 'Incomplete MANE identity')
            _require(re.fullmatch(r'ENST\d+\.\d+',r['transcript_id']) is not None, 'Versioned MANE transcript required')
            _require(r['transcript_id'] not in seen and r['strand'] in ('+','-') and bool(r['gene_id']) and bool(r['gene_symbol']), 'Ambiguous MANE identity')
            seen.add(r['transcript_id'])
        object.__setattr__(self, 'data', _frozen(d))

    @classmethod
    def from_bytes(cls, payload, *, expected_sha256):
        require_sha(expected_sha256)
        _require(len(payload) <= 16 * 1024**2 and digest(payload) == expected_sha256, 'Source manifest SHA256 mismatch or oversized manifest')
        return cls(json.loads(payload), expected_sha256)


@dataclass(frozen=True)
class SpliceScore:
    gene_symbol: str
    scores: tuple
    offsets: tuple


@dataclass(frozen=True)
class SourceResult:
    records: tuple[SpliceScore, ...]
    artifact_bytes: bytes
    receipts: tuple = ()


def parse_vcf_scores(payload, *, chromosome, position_1based, ref, alt):
    """Select exact genomic allele; dots stay missing and zeros stay numeric."""
    _require(type(position_1based) is int and position_1based > 0, 'Positive VCF coordinate required')
    choices = set()
    for line in payload.decode('utf-8').splitlines():
        if not line or line.startswith('#'): continue
        cells = line.split('\t')
        _require(len(cells) >= 8, 'Malformed VCF row')
        if cells[0].removeprefix('chr') != chromosome.removeprefix('chr') or int(cells[1]) != position_1based: continue
        _require(cells[3] == ref, 'Source REF conflicts with normalized reference allele')
        if alt not in cells[4].split(','): continue
        info = [v.partition('=')[2] for v in cells[7].split(';') if v.startswith('SpliceAI=')]
        _require(len(info) <= 1, 'Duplicate SpliceAI field')
        for annotation in info[0].split(',') if info and info[0] != '.' else ():
            f = annotation.split('|')
            _require(len(f) == 10, 'Unexpected SpliceAI schema')
            _require(f[0] in cells[4].split(','), 'Annotation allele absent from VCF ALT')
            if f[0] != alt: continue
            scores = tuple(None if v == '.' else float(v) for v in f[2:6])
            offsets = tuple(None if v == '.' else int(v) for v in f[6:10])
            _require(bool(f[1]) and all(v is None or math.isfinite(v) and 0 <= v <= 1 for v in scores), 'Invalid SpliceAI score')
            choices.add(SpliceScore(f[1],scores,offsets))
    return tuple(sorted(choices,key=repr))


def tabix_chunks(index_bytes, chromosome, position_1based, *, total_bytes):
    """Read a gzip TBI VCF index, selecting standard reg2bins at one base."""
    try:
        data = _gunzip(index_bytes,64 * 1024**2)
        _require(data[:4] == b'TBI\1', 'Invalid TBI magic')
        n,fmt,seq,beg,end,meta,skip,names_len = struct.unpack_from('<8i',data,4)
        _require((fmt,seq,beg,end,meta,skip) == (2,1,2,0,35,0) and 0 < n < 10000 and 0 < names_len <= len(data)-36,'Unexpected VCF index schema')
        offset = 36 + names_len
        names = data[36:offset].rstrip(b'\0').decode().split('\0')
        _require(len(names) == n and len(set(names)) == n,'Invalid TBI names')
        zero=position_1based-1
        _require(type(position_1based) is int and 0 <= zero < 2**29,'Position outside TBI support')
        bins_wanted = (0,1+(zero>>26),9+(zero>>23),73+(zero>>20),585+(zero>>17),4681+(zero>>14))
        spans=[]
        for name in names:
            count,=struct.unpack_from('<i',data,offset); offset+=4
            _require(0 <= count <= 37451,'Invalid TBI bin count')
            bins={}
            for _ in range(count):
                b,k=struct.unpack_from('<Ii',data,offset); offset+=8
                _require(0 <= k <= (len(data)-offset)//16 and b not in bins,'Invalid TBI chunks')
                bins[b]=tuple(struct.iter_unpack('<QQ',data[offset:offset+16*k])); offset+=16*k
            k,=struct.unpack_from('<i',data,offset); offset+=4
            _require(0 <= k <= (len(data)-offset)//8,'Invalid TBI linear index')
            if name.removeprefix('chr') == chromosome.removeprefix('chr'):
                linear=zero>>14
                minimum=struct.unpack_from('<Q',data,offset+8*linear)[0] if linear < k else 0
                for b in bins_wanted:
                    for lo,hi in bins.get(b,()):
                        _require(lo < hi and hi>>16 < total_bytes,'Invalid TBI virtual offset')
                        if hi > minimum: spans.append((max(lo,minimum),hi))
            offset+=8*k
        _require(len(data)-offset in (0,8),'Unexpected TBI trailer')
        merged=[]
        for lo,hi in sorted(set(spans)):
            if merged and lo>>16 <= merged[-1][1]>>16: merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
            else: merged.append((lo,hi))
        return tuple(merged)
    except (struct.error,UnicodeError,EOFError,zlib.error,OSError) as exc:
        raise ValueError('Invalid or truncated TBI') from exc


def range_for_chunk(span, *, total_bytes):
    lo,hi=span
    _require(0 <= lo < hi and hi>>16 < total_bytes,'Invalid virtual span')
    start=lo>>16
    end=min(total_bytes,(hi>>16)+(65536 if hi & 65535 else 0))
    _require(start < end,'Empty compressed range')
    return start,end-1


def selected_bgzf(body,span):
    """Validate block CRC/ISIZE and return exactly one tabix virtual span."""
    lo,hi=span; base=lo>>16; at=0; pieces=[]
    try:
        while at < len(body) and base+at <= hi>>16:
            address=base+at
            if address == hi>>16 and not hi & 65535: break
            _require(body[at:at+4] == b'\x1f\x8b\x08\x04' and body[at+12:at+16] == b'BC\x02\x00','Invalid BGZF header')
            _require(struct.unpack_from('<H',body,at+10)[0] == 6,'Unqualified BGZF extra fields')
            length=struct.unpack_from('<H',body,at+16)[0]+1
            _require(28 <= length <= 65536 and at+length <= len(body),'Truncated BGZF block')
            # Gzip ISIZE supplies a bounded precondition before decompression.
            size=struct.unpack_from('<I',body,at+length-4)[0]
            _require(size <= 65536,'Oversized BGZF uncompressed block')
            decoder=zlib.decompressobj(31)
            decoded=decoder.decompress(body[at:at+length],65537)
            _require(len(decoded) == size and decoder.eof and not decoder.unconsumed_tail,'Invalid BGZF decoded length')
            start=lo & 65535 if at == 0 else 0
            end=hi & 65535 if address == hi>>16 else len(decoded)
            _require(start <= end <= len(decoded),'Invalid BGZF virtual slice')
            pieces.append(decoded[start:end]); at+=length
        _require(base+at >= hi>>16,'Incomplete BGZF span')
        return b''.join(pieces)
    except (struct.error,zlib.error) as exc: raise ValueError('Invalid BGZF block') from exc


class LocalSpliceSource:
    def __init__(self,manifest,path,*,max_bytes=16*1024**2):
        self.manifest=manifest; self.path=Path(path); self.max_bytes=max_bytes
        require_sha(manifest.data.get('score_sha256'))

    def query(self,chromosome,position_1based,ref,alt,**kwargs):
        with self.path.open('rb') as stream: raw=stream.read(self.max_bytes+1)
        _require(len(raw) <= self.max_bytes and digest(raw) == self.manifest.data['score_sha256'],'Score source SHA256 mismatch or oversized source')
        payload=_gunzip(raw,self.max_bytes) if raw.startswith(b'\x1f\x8b') else raw
        records=parse_vcf_scores(payload,chromosome=chromosome,position_1based=position_1based,ref=ref,alt=alt)
        return SourceResult(records,raw,({'source_sha256':digest(raw)},))


class IndexedSpliceSource:
    def __init__(self,manifest,fetcher,*,max_index_bytes=8*1024**2,max_range_bytes=16*1024**2):
        self.manifest=manifest; self.fetcher=fetcher
        self.max_index_bytes=max_index_bytes; self.max_range_bytes=max_range_bytes
        _require(all(k in manifest.data for k in MANIFEST_SCHEMA['indexed']),'Incomplete indexed source manifest')
        self._index=None

    def query(self,chromosome,position_1based,ref,alt,*,variant_key='spliceai'):
        d=self.manifest.data
        if self._index is None:
            raw=self.fetcher.fetch(d['index_url'],variant_key=variant_key,expected_etag=d.get('index_etag'))
            _require(len(raw) <= self.max_index_bytes and digest(raw) == d['index_sha256'],'TBI SHA256 mismatch or oversized index')
            self._index=raw
        spans=tabix_chunks(self._index,chromosome,position_1based,total_bytes=d['score_bytes'])
        ranges=[range_for_chunk(s,total_bytes=d['score_bytes']) for s in spans]
        _require(sum(b-a+1 for a,b in ranges) <= self.max_range_bytes,'Indexed request exceeds source range bound')
        pieces=[]; receipts=[]
        for span,(start,end) in zip(spans,ranges):
            body=self.fetcher.fetch(d['score_url'],variant_key=variant_key,expected_etag=d['score_etag'],byte_range=(start,end))
            _require(len(body) == end-start+1,'Range response has wrong length; whole-file fallback prohibited')
            piece=selected_bgzf(body,span)
            _require(not piece or piece.endswith(b'\n'),'TBI span does not end on VCF row boundary')
            pieces.append(piece); receipts.append({'start':start,'end_inclusive':end,'sha256':digest(body)})
        payload=b''.join(pieces)
        records=parse_vcf_scores(payload,chromosome=chromosome,position_1based=position_1based,ref=ref,alt=alt)
        return SourceResult(records,payload,tuple(receipts))
