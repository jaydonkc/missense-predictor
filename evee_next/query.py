"""Build label-free one-case inference inputs from a pinned missense mapping."""

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import sqlite3
from typing import Mapping

from .reference import ReferenceResolver, STUDY, load_verified_module
from .report import MissenseMapping, NormalizedVariant, VariantRequest, matching_mappings


QUERY_SCOPE = "evee_next_query_missense_v1"
WINDOW_CONTRACT_SHA256 = "9077abd19e5e387f57ec07470fa1c5ef62143224a748330b443c70ce4468d41a"
WINDOW_CONTRACT = Path(__file__).resolve().parents[1] / STUDY / "model_intake" / "window_contract.py"
DNA_WIDTH = 8192
PROTEIN_MAXIMUM = 1022
MINIMUM_DIRECTIONAL_CONTEXT = 2049


@dataclass(frozen=True)
class PreparedQuery:
    variant: NormalizedVariant
    mapping: MissenseMapping
    payload: bytes
    receipt: dict

    def write(self, directory: Path) -> None:
        """Publish one query in a new directory, writing its receipt last."""
        directory.mkdir()
        with (directory / "query_paired_inputs.json").open("xb") as stream:
            stream.write(self.payload)
        with (directory / "preparation_receipt.json").open("x", encoding="utf-8") as stream:
            json.dump(self.receipt, stream, indent=2, allow_nan=False)
            stream.write("\n")


@dataclass(frozen=True)
class UnpreparedQuery:
    status: str
    detail: str
    variant: NormalizedVariant | None
    mapping_options: tuple[MissenseMapping, ...]
    provenance: Mapping[str, str]


def sequence_sha(sequence):
    return sha256(sequence.encode("ascii")).hexdigest()


def prepare_query(request: VariantRequest, resolver: ReferenceResolver) -> PreparedQuery | UnpreparedQuery:
    """Prepare both input modalities only after a unique exact transcript/protein match."""
    resolution = resolver.resolve(request)

    def failed(status, detail, options=resolution.mappings):
        return UnpreparedQuery(status, detail, resolution.variant, options, resolution.provenance)

    variant, options = resolution.variant, resolution.mappings
    if variant is None or not options:
        return failed(resolution.status, resolution.detail)
    matched = matching_mappings(request, options)
    if not matched:
        return failed("mapping_conflict", "The requested transcript/protein versions do not match any exact mapping.")
    if len(matched) != 1:
        return failed("ambiguous_mapping", "Choose a versioned transcript or protein from the exact mapping options.", matched)
    mapping, = matched
    try:
        contract = load_verified_module("_evee_query_window_contract", WINDOW_CONTRACT, WINDOW_CONTRACT_SHA256)
        with resolver._query_sources(variant, mapping) as (genome, protein):
            dna = contract.centered_window(genome.length, variant.position_0based, DNA_WIDTH)
            length = dna.end_0based_exclusive - dna.start_0based
            index = dna.variant_index_0based
            if length != DNA_WIDTH or min(index, length - index - 1) < MINIMUM_DIRECTIONAL_CONTEXT:
                return failed("outside_frozen_feature_scope", "Both DNA orientations need a complete 8,192-base window and "
                              "nonempty distal causal pools.", matched)
            genome_ref = genome.read(dna.start_0based, dna.end_0based_exclusive)
            # DNA stays in forward genomic orientation. The resolver already checked strand-specific codon translation.
            dna_pair = contract.dna_pair(genome_ref, index, variant.ref, variant.alt)
            protein_pair = contract.protein_pair(protein, mapping.residue_1based, mapping.ref_aa, mapping.alt_aa,
                                                 maximum_length=PROTEIN_MAXIMUM)
    except FileNotFoundError as error:
        return failed("missing_reference", str(error), matched)
    except (OSError, ValueError, sqlite3.Error) as error:
        return failed("reference_integrity_error", str(error), matched)

    case = {"variant_id": variant.variant_id, "protein_id": mapping.protein_id, "transcript_id": mapping.transcript_id,
            "protein_ref": protein_pair["ref_window"], "protein_alt": protein_pair["alt_window"],
            "protein_variant_index0": protein_pair["variant_index_0based"],
            "dna_windows": [{"length": length, "variant_index0": index, "ref": dna_pair["fwd_ref"], "alt": dna_pair["fwd_alt"]}]}
    payload = (json.dumps({"scope": QUERY_SCOPE, "cases": [case]}, separators=(",", ":"), sort_keys=True) + "\n").encode("utf-8")
    receipt = {"schema": "evee-next-query-preparation/v1", "status": "prepared_label_free_query", "scope": QUERY_SCOPE,
               "payload_sha256": sha256(payload).hexdigest(), "variant": asdict(variant), "mapping": asdict(mapping),
               "reference_provenance": dict(resolution.provenance), "window_contract_sha256": WINDOW_CONTRACT_SHA256,
               "dna_window": {**asdict(dna), "length": length,
                              "reverse_complement_variant_index_0based": dna_pair["revcomp_variant_index_0based"],
                              "forward_ref_sha256": sequence_sha(dna_pair["fwd_ref"]),
                              "forward_alt_sha256": sequence_sha(dna_pair["fwd_alt"]),
                              "reverse_complement_ref_sha256": sequence_sha(dna_pair["revcomp_ref"]),
                              "reverse_complement_alt_sha256": sequence_sha(dna_pair["revcomp_alt"])},
               "protein_window": {key: protein_pair[key] for key in ("start_0based", "end_0based_exclusive",
                                  "variant_index_0based", "esm_token_index_0based")} | {
                   "source_sequence_sha256": sequence_sha(protein), "ref_sha256": sequence_sha(protein_pair["ref_window"]),
                   "alt_sha256": sequence_sha(protein_pair["alt_window"])}}
    return PreparedQuery(variant, mapping, payload, receipt)
