"""Attach exact protein-matched cached effect diagnostics to genomic reports.

This adapter does not change a clinical score. The caller supplies the existing
standalone VariantAnalyzer/ReferenceResolver, retaining nucleotide/codon checks.
"""
from .pipeline import ContextRecord, EffectRecord
from .runtime import CachedEffectProvider, PROTEINS


class GenomicCachedBiology:
    def __init__(self,provider: CachedEffectProvider):
        self.provider=provider
        self.genes={protein:gene for gene,protein in PROTEINS.items()}

    def _report(self,mapping):
        gene=self.genes.get(mapping.protein_id)
        if gene is None:return None
        return self.provider.query(gene,mapping.residue_1based,mapping.ref_aa,mapping.alt_aa)['molecular_effect']

    def context(self,mapping):
        r=self._report(mapping)
        if r is None:return ContextRecord('unavailable')
        ctx=r['context']
        values={k:ctx[k] for k in ('reference_pLDDT','reference_RSA')}
        return ContextRecord('available' if all(v is not None for v in values.values()) else 'partial',values,{'context':ctx['source'],'result':r['source']})

    def effects(self,mapping):
        r=self._report(mapping)
        if r is None:return ()
        return (EffectRecord(r['endpoint'],r['predictions']['J'],r['units'],r['source'],r['qualification'],r['uncertainty']),)
