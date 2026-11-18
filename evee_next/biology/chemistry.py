"""Directed substitution properties from the project's pinned classical tables."""
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path

AA = 'ACDEFGHIKLMNPQRSTVWY'
TABLE_PATH = Path('research/missense-classical-tree-2026-10-26/sources/tables.json')
TABLE_SHA = '134751dbc64825f0614b9490fdc98f67b40f9da621dffd9d8c64d7eee1cd768a'
PROPERTIES = ('charge', 'volume', 'hydrophobicity', 'polarity')
NUMERIC_COLUMNS = ('grantham', 'pam250') + tuple('delta_'+x for x in PROPERTIES) + tuple('abs_delta_'+x for x in PROPERTIES) + tuple(x+'_'+a for a in 'GPC' for x in ('lose','gain'))


@dataclass(frozen=True)
class ProteinVariant:
    gene: str
    protein_id: str
    position_1based: int
    ref_aa: str
    alt_aa: str
    reference_sha256: str

    def validate(self, sequence: str) -> None:
        if (not self.gene or not self.protein_id or type(self.position_1based) is not int
                or not 1 <= self.position_1based <= len(sequence)
                or len(self.ref_aa)!=1 or self.ref_aa not in AA
                or len(self.alt_aa)!=1 or self.alt_aa not in AA or self.ref_aa==self.alt_aa
                or sequence[self.position_1based-1]!=self.ref_aa
                or hashlib.sha256(sequence.encode()).hexdigest()!=self.reference_sha256):
            raise ValueError('Protein identity, reference sequence or substitution conflict')


class Chemistry:
    def __init__(self, tables: dict, source: str):
        self.tables = tables
        self.source = source

    @classmethod
    def from_project(cls, root: Path):
        payload=(root/TABLE_PATH).read_bytes()
        if hashlib.sha256(payload).hexdigest()!=TABLE_SHA:
            raise ValueError('Chemistry table hash differs')
        return cls(json.loads(payload), str(TABLE_PATH)+'#sha256='+TABLE_SHA)

    def calculate(self, ref: str, alt: str) -> dict:
        if len(ref)!=1 or ref not in AA or len(alt)!=1 or alt not in AA or ref==alt:
            raise ValueError('Canonical missense substitution required')
        values={k:float(self.tables[k][ref][alt]) for k in ('grantham','pam250')}
        for name in PROPERTIES:
            delta=float(self.tables[name][alt])-float(self.tables[name][ref])
            values['delta_'+name]=delta
            values['abs_delta_'+name]=abs(delta)
        for a in 'GPC':
            values['lose_'+a]=float(ref==a)
            values['gain_'+a]=float(alt==a)
        if not all(math.isfinite(x) for x in values.values()):
            raise ValueError('Nonfinite chemical property')
        return dict(values, source=self.source, is_pathogenicity=False,
                    interpretation='chemical_difference_not_functional_damage')

    def vector(self, ref: str, alt: str):
        values=self.calculate(ref,alt)
        return [values[k] for k in NUMERIC_COLUMNS]+[float(ref==a) for a in AA]+[float(alt==a) for a in AA]
