"""Attach the two reviewed RNA observations independently of model scoring."""

from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path

from .report import EvidenceRecord, VariantEvidenceReport


LEDGER = Path(
    "research/missense-benchmark-2026-10-23/clinical_assertion_review/scaleup_v2/"
    "primary_evidence_v2/incremental_packet_v7/batch/RNA_product_mismatch_ledger.json"
)
LEDGER_SHA256 = "6e5ebeaaa2a7d44d32a744fcd88e613706dfc88cef3f6120d5702c972b7be845"

CASE_CONTEXT = {
    "chr10:123053159:G:A": (
        "The study observed exon 10 skipping in fibroblast cDNA from a homozygous proband. "
        "Its biochemical experiment assayed the splice product. The fraction of normally spliced RNA "
        "and neurological penetrance remain unresolved."
    ),
    "chr3:10090314:G:A": (
        "The study observed use of a cryptic splice acceptor in a patient who also carried c.2835dupC. "
        "The fraction of aberrant RNA and the effect of the isolated amino-acid substitution remain unresolved."
    ),
}


def attach_reviewed_rna_evidence(report: VariantEvidenceReport, project_root: Path) -> VariantEvidenceReport:
    """Add case evidence for an exact GRCh38 allele in a resolved report."""
    variant = report.normalized_variant
    if (report.mapping is None or variant is None or variant.assembly != "GRCh38"
            or variant.organism != "Homo sapiens" or variant.variant_id not in CASE_CONTEXT):
        return report

    ledger = project_root / LEDGER
    try:
        payload = ledger.read_bytes()
    except FileNotFoundError:
        return report
    if sha256(payload).hexdigest() != LEDGER_SHA256:
        raise ValueError(f"Reviewed RNA ledger SHA-256 mismatch: {ledger}")

    case = next(record for record in json.loads(payload)["records"] if record["variant_id"] == variant.variant_id)
    start, end = case["reported_RNA_deletion_coding_coordinates"]
    observation = EvidenceRecord(
        "reported_effect",
        f"Published RNA observation for {case['gene']} {case['coding_HGVS']}: "
        f"{case['reported_RNA_deleted_nt']}-nt deletion spanning coding positions {start}–{end} "
        f"on this RefSeq transcript. {CASE_CONTEXT[variant.variant_id]} "
        "This case evidence is separate from the scoring model.",
        case["source_proofs"][0]["url"],
    )
    if observation in report.evidence:
        return report
    return replace(report, evidence=(*report.evidence, observation), provenance={
        **report.provenance,
        "reviewed_rna_ledger": str(ledger.resolve()),
        "reviewed_rna_ledger_sha256": LEDGER_SHA256,
    })
