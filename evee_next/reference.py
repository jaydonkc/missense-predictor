"""Offline Ensembl 116 missense identity resolution against pinned GRCh38 sources."""

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from contextlib import closing, contextmanager
import csv
from dataclasses import dataclass
from hashlib import file_digest, sha256
import importlib.util
from io import StringIO
from itertools import accumulate
import json
import os
from pathlib import Path
import re
import sqlite3
from threading import RLock
from types import MappingProxyType
from typing import Mapping

from .report import AA, CHROMOSOMES, IdentityResolution, MissenseMapping, NormalizedVariant, VariantRequest, normalize


STUDY = Path("research/missense-benchmark-2026-10-23")
CONCORDANCE_RECEIPT_SHA256 = "92a2c37f05bbb9cd71db25904107a66c7b165fc93fad2641c46cb373e4dad162"
GENOME_MANIFEST_SHA256 = "2f929a0487f87c59329a80e9aa6c3a80688de14f0f02347f7f7b31e898e1fc7c"
TRANSLATION_SOURCE_SHA256 = "3ab1008d8f15f796a326031f09402783de7bf93e2651e3460f8ac4ff54f76a6c"
EXACT = "exact_full_cds_genome_match"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
CONSEQUENCES = {"initiation_context_requires_review", "synonymous", "stop_gained", "stop_lost", "stop_retained"}


def check_hash(path, expected):
    with path.open("rb") as stream:
        actual = file_digest(stream, "sha256").hexdigest()
    if actual != expected:
        raise ValueError(f"Pinned reference SHA-256 mismatch: {path.name}")


def verified_bytes(path, expected):
    data = path.read_bytes()
    if sha256(data).hexdigest() != expected:
        raise ValueError(f"Pinned reference SHA-256 mismatch: {path.name}")
    return data


def load_verified_module(name, path, expected):
    source = verified_bytes(path, expected)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    exec(compile(source, str(path), "exec"), vars(module))
    return module


def signature(metadata):
    return metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns


@dataclass(frozen=True)
class TranscriptAssessment:
    transcript_id: str
    protein_id: str | None
    status: str
    detail: str
    reference_codon: str | None = None
    alternate_codon: str | None = None
    cds_offset_0based: int | None = None
    mapping: MissenseMapping | None = None


@dataclass(frozen=True)
class TranscriptContext:
    """Qualified reference context; CDS-only intervals are never exon boundaries."""
    transcript_id: str
    protein_id: str | None
    strand: str
    exon_intervals: tuple[tuple[int, int], ...]
    cds_intervals: tuple[tuple[int, int], ...]
    reference_hashes: Mapping[str, str]
    boundary_status: str
    protein_sequence: str | None = None

    def __post_init__(self):
        if self.strand not in {'+','-'} or self.boundary_status not in {'complete_exons','cds_only','unavailable'}:
            raise ValueError('Explicit transcript strand and boundary completeness required')
        for name in ('exon_intervals','cds_intervals'):
            rows=tuple(tuple(row) for row in getattr(self,name))
            if any(len(row)!=2 or any(type(v) is not int for v in row) or not 0<=row[0]<row[1] for row in rows):
                raise ValueError('Intervals must be zero-based half-open')
            if rows != tuple(sorted(rows)) or any(a[1]>b[0] for a,b in zip(rows,rows[1:])):
                raise ValueError('Reference intervals must be ordered and nonoverlapping')
            object.__setattr__(self,name,rows)
        if self.boundary_status=='complete_exons' and not self.exon_intervals:
            raise ValueError('Complete exon context requires qualified exons')
        object.__setattr__(self,'reference_hashes',MappingProxyType(dict(self.reference_hashes)))


@dataclass(frozen=True)
class ReferenceResolution:
    status: str
    detail: str
    provenance: Mapping[str, str]
    variant: NormalizedVariant | None = None
    assessments: tuple[TranscriptAssessment, ...] = ()

    @property
    def mappings(self) -> tuple[MissenseMapping, ...]:
        return tuple(row.mapping for row in self.assessments if row.mapping is not None)


class ReferenceResolver:
    """Resolve within the selected transcript index; caller supplies immutable FASTA/FAI files."""

    @classmethod
    def from_project(cls, root: Path, fasta_directory: Path):
        study = root.resolve() / STUDY
        return cls(study / "reference_bundle", study / "genome_reference", fasta_directory,
                   concordance_receipt_sha256=CONCORDANCE_RECEIPT_SHA256, manifest_sha256=GENOME_MANIFEST_SHA256,
                   translation_sha256=TRANSLATION_SOURCE_SHA256)

    def __init__(self, bundle: Path, genome_sources: Path, fasta_directory: Path, *, concordance_receipt_sha256: str,
                 manifest_sha256: str, translation_sha256: str):
        receipt_path = bundle / "genome_validation" / "cds_genome_receipt.json"
        manifest_path = genome_sources / "remote_manifest.json"
        receipt = json.loads(verified_bytes(receipt_path, concordance_receipt_sha256))
        manifest = json.loads(verified_bytes(manifest_path, manifest_sha256))
        if (receipt["status"] != "completed_full_cds_genome_audit" or receipt["reference_release"] != 116
                or receipt["reference_manifest_sha256"] != manifest_sha256 or receipt["outcomes_read"] is not False
                or manifest["status"] != "complete_canonical_reference" or manifest["release"] != 116
                or manifest["assembly"] != "GRCh38" or manifest["canonical_nuclear_chromosomes_only"] is not True):
            raise ValueError("Expected the pinned Ensembl 116 GRCh38 full CDS concordance bundle")
        self._index = (bundle / "selected_reference_index.sqlite").resolve(strict=True)
        self._index_state = signature(self._index.stat())
        check_hash(self._index, receipt["index_sha256"])
        if signature(self._index.stat()) != self._index_state:
            raise ValueError("Reference index changed during verification")
        table = receipt_path.with_name("cds_genome_concordance.tsv")
        table_bytes = verified_bytes(table, receipt["output_sha256"])
        audit_rows = list(csv.DictReader(StringIO(table_bytes.decode("utf-8"), newline=""), delimiter="\t"))
        self._audit = {(row["transcript_id"], row["protein_id"]): row for row in audit_rows}
        if len(self._audit) != len(audit_rows) or len(audit_rows) != receipt["rows"] or Counter(
                row["status"] for row in audit_rows) != receipt["status_counts"]:
            raise ValueError("Concordance table disagrees with the full CDS audit receipt")
        reader = load_verified_module("_evee_indexed_reference", genome_sources / "verify_clinical_ref.py", receipt["reader_sha256"])
        merge = load_verified_module("_evee_full_cds_audit", bundle / "verify_cds_against_genome.py", receipt["code_sha256"])
        translation = load_verified_module("_evee_translation", bundle / "translation_and_codon_audit.py", translation_sha256)
        self._indexed_chromosome = reader.IndexedChromosome
        self._merge = merge.merged_intervals
        self._translation = translation.translation
        self._genetic_code = translation.STANDARD
        self._fasta_directory = fasta_directory.resolve()
        self._sources = {row["chromosome"]: row["fasta"] for row in manifest["files"]}
        if len(self._sources) != len(manifest["files"]) or not set(self._sources) <= CHROMOSOMES:
            raise ValueError("Duplicate or noncanonical chromosome entries in the reference manifest")
        self._chromosomes = {}
        self._lock = RLock()
        self._closed = False
        self._protein_ids = defaultdict(list)
        indexed = defaultdict(list)
        with closing(self._connect()) as db:
            for tid, pid in db.execute("SELECT transcript_id,protein_id FROM protein ORDER BY protein_id"):
                self._protein_ids[tid].append(pid)
            for chromosome, start, end, tid in db.execute("""
                SELECT chromosome,start_1based,end_1based,transcript_id FROM feature
                WHERE feature IN ('CDS','start_codon','stop_codon')
                ORDER BY chromosome,start_1based,end_1based,transcript_id
            """):
                indexed[chromosome].append((start, end, tid))
        self._regions = {chrom: (tuple(x[0] for x in rows), tuple(accumulate((x[1] for x in rows), max)), tuple(rows))
                         for chrom, rows in indexed.items()}
        self.provenance = {"identity_mapping": str(receipt_path.resolve()), "reference_release": "Ensembl 116 / GRCh38.p14",
                           "identity_scope": "selected_ensembl_transcript_index", "reference_index_sha256": receipt["index_sha256"],
                           "cds_genome_receipt_sha256": concordance_receipt_sha256,
                           "cds_genome_concordance_sha256": receipt["output_sha256"], "genome_manifest_sha256": manifest_sha256,
                           "indexed_reader_sha256": receipt["reader_sha256"], "concordance_method_sha256": receipt["code_sha256"],
                           "translation_source_sha256": translation_sha256}

    def _connect(self):
        if signature(self._index.stat()) != self._index_state:
            raise ValueError("Pinned reference index changed after verification")
        return sqlite3.connect(self._index.as_uri() + "?mode=ro&immutable=1", uri=True)

    def _chromosome(self, name):
        source = self._sources[name]
        fasta_name, fai_name = source["path"], source["index"]
        if Path(fasta_name).name != fasta_name or Path(fai_name).name != fai_name:
            raise ValueError("Reference paths must be files in the supplied FASTA directory")
        fasta, fai = (self._fasta_directory / fasta_name, self._fasta_directory / fai_name)
        states = signature(fasta.stat()), signature(fai.stat())
        if name not in self._chromosomes:
            check_hash(fasta, source["sha256"])
            check_hash(fai, source["index_sha256"])
            reader = self._indexed_chromosome(fasta, fai)
            if (reader.name != name or reader.length != source["length"] or signature(os.fstat(reader.file.fileno())) != states[0]
                    or (signature(fasta.stat()), signature(fai.stat())) != states):
                reader.close()
                raise ValueError("Indexed chromosome metadata or reference bytes changed")
            self._chromosomes[name] = reader, states
        reader, original_states = self._chromosomes[name]
        if states != original_states:
            raise ValueError("Pinned chromosome source changed after verification")
        return reader

    def _candidate_transcripts(self, chromosome, position):
        starts, prefix_ends, regions = self._regions.get(chromosome, ((), (), ()))
        start = bisect_left(prefix_ends, position)
        end = bisect_right(starts, position)
        return sorted({tid for _, stop, tid in regions[start:end] if stop >= position})

    def _assess(self, db, genome, variant, tid, pid):
        def result(status, detail, *extra):
            return TranscriptAssessment(tid, pid, status, detail, *extra)

        if (not re.fullmatch(r"ENST\d+\.\d+", tid or "") or not re.fullmatch(r"ENSP\d+\.\d+", pid or "")):
            return result("missing_versioned_identity", "No exact versioned ENST/ENSP pair is available.")
        evidence = self._audit.get((tid, pid))
        if evidence is None or evidence["status"] != EXACT:
            return result("full_cds_concordance_unavailable", "Full CDS-genome concordance is absent or failed.")
        row = db.execute("SELECT p.sequence,c.sequence FROM protein p JOIN cds c USING(transcript_id) "
                         "WHERE p.transcript_id=? AND p.protein_id=?", (tid, pid)).fetchone()
        if row is None:
            return result("missing_cds_or_protein", "The pinned index lacks a CDS or protein sequence for this mapping.")
        protein, cds = row
        features = db.execute("SELECT chromosome,strand,feature,start_1based,end_1based,phase,protein_id FROM feature "
                              "WHERE transcript_id=? AND feature IN ('CDS','stop_codon','start_codon')", (tid,)).fetchall()
        strand = evidence["strand"]
        location = {(row[0], row[1]) for row in features}
        if (strand not in {"+", "-"} or evidence["chromosome"] != variant.chromosome
                or location != {(variant.chromosome, strand)} or any(not 1 <= s <= e <= genome.length for _, _, _, s, e, _, _ in features)):
            return result("ambiguous_feature_location", "CDS annotation and concordance disagree on chromosome, strand or bounds.")
        coding = sorted(((s, e, p, protein_id) for _, _, f, s, e, p, protein_id in features if f == "CDS"),
                        reverse=strand == "-")
        cumulative = 0
        if not coding:
            return result("missing_cds_annotation", "No CDS segments are available for this versioned transcript.")
        for start, end, phase, protein_id in coding:
            if phase != str((-cumulative) % 3) or protein_id != pid:
                return result("unsupported_cds_annotation", "CDS phases or protein annotations do not support an exact complete mapping.")
            cumulative += end - start + 1
        ascending = sorted(coding)
        if any(left[1] >= right[0] for left, right in zip(ascending, ascending[1:])):
            return result("unsupported_cds_annotation", "Overlapping CDS segments give ambiguous coordinates.")

        spans = self._merge([(s, e) for _, _, f, s, e, _, _ in features if f in {"CDS", "stop_codon"}])
        forward = "".join(genome.read(start - 1, end) for start, end in spans)
        reconstructed = forward if strand == "+" else forward.translate(COMPLEMENT)[::-1]
        if (reconstructed != cds or evidence["mismatch_bases"] != "0" or int(evidence["interval_count"]) != len(spans)
                or int(evidence["cds_reference_length"]) != len(cds) or int(evidence["genome_reconstructed_length"]) != len(cds)
                or sha256(cds.encode()).hexdigest() != evidence["cds_sha256"]
                or evidence["cds_sha256"] != evidence["genome_reconstructed_sha256"]):
            return result("cds_genome_disagreement", "The indexed genomic intervals disagree with the pinned full CDS audit.")
        translation_status, _ = self._translation(cds, protein)
        if not protein or not set(protein) <= set(AA) or translation_status not in {
                "exact_translation", "exact_with_documented_initiator_methionine"}:
            return result("unsupported_translation", f"The CDS lacks an exact supported protein translation ({translation_status}).")

        pos = variant.position_1based
        offset = 0
        for start, end in (spans if strand == "+" else reversed(spans)):
            if start <= pos <= end:
                offset += pos - start if strand == "+" else end - pos
                break
            offset += end - start + 1
        else:
            return result("outside_cds", "The queried site is outside the reconstructed CDS and terminal stop.")
        ref, alt = (variant.ref, variant.alt) if strand == "+" else (variant.ref.translate(COMPLEMENT), variant.alt.translate(COMPLEMENT))
        start = offset - offset % 3
        reference_codon = cds[start:start + 3]
        i = offset % 3
        alternate_codon = reference_codon[:i] + alt + reference_codon[i + 1:]
        codons = (reference_codon, alternate_codon, offset)
        if cds[offset] != ref:
            return result("cds_reference_mismatch", "The strand-adjusted genomic REF differs from the transcript CDS.", *codons)
        ref_aa, alt_aa = self._genetic_code[reference_codon], self._genetic_code[alternate_codon]
        residue = offset // 3 + 1
        if residue == 1:
            return result("initiation_context_requires_review", "The substitution affects the first codon; initiation requires separate review.", *codons)
        if ref_aa == "*":
            return result("stop_retained" if alt_aa == "*" else "stop_lost", "The substitution affects the terminal stop codon.", *codons)
        if alt_aa == "*":
            return result("stop_gained", "The substitution introduces a stop codon.", *codons)
        if ref_aa == alt_aa:
            return result("synonymous", "The substitution preserves the translated amino acid.", *codons)
        mapping = MissenseMapping(tid, pid, strand, residue, ref_aa, alt_aa)
        return result("exact_missense", "The pinned transcript, protein and genomic sequence agree.", *codons, mapping)

    def resolve(self, value: VariantRequest | NormalizedVariant) -> ReferenceResolution:
        provenance = dict(self.provenance)
        try:
            if isinstance(value, NormalizedVariant):
                request = VariantRequest(value.chromosome, value.position_1based, value.ref, value.alt,
                                         organism=value.organism, assembly=value.assembly)
                variant = normalize(request)
                if variant != value:
                    raise ValueError("Inconsistent normalized variant coordinates or allele identity")
            else:
                variant = normalize(value)
        except NotImplementedError as error:
            return ReferenceResolution("unsupported", str(error), provenance)
        except ValueError as error:
            return ReferenceResolution("invalid_request", str(error), provenance)

        with self._lock:
            if self._closed:
                raise RuntimeError("Reference resolver is closed")
            if variant.chromosome not in self._sources:
                return ReferenceResolution("missing_reference", "No pinned FASTA and index are available for this chromosome.",
                                           provenance, variant)
            try:
                genome = self._chromosome(variant.chromosome)
                source = self._sources[variant.chromosome]
                provenance.update(chromosome_fasta_sha256=source["sha256"], chromosome_fai_sha256=source["index_sha256"])
                if variant.position_1based > genome.length:
                    return ReferenceResolution("invalid_request", "Position is beyond the pinned chromosome length.", provenance, variant)
                base = genome.read(variant.position_0based, variant.position_1based)
                if base == "N":
                    return ReferenceResolution("missing_reference", "The pinned genome has an ambiguous base at this position.", provenance, variant)
                if base != variant.ref:
                    return ReferenceResolution("reference_mismatch", f"The pinned genomic REF at this position is {base}.", provenance, variant)
                candidates = self._candidate_transcripts(variant.chromosome, variant.position_1based)
                with closing(self._connect()) as db:
                    assessments = tuple(self._assess(db, genome, variant, tid, pid) for tid in candidates
                                        for pid in self._protein_ids.get(tid, (None,)))
            except FileNotFoundError as error:
                return ReferenceResolution("missing_reference", f"Pinned indexed reference is unavailable: {error}", provenance, variant)
            except (KeyError, OSError, ValueError, sqlite3.Error) as error:
                return ReferenceResolution("reference_integrity_error", str(error), provenance, variant)
        exact = sum(row.mapping is not None for row in assessments)
        if exact:
            return ReferenceResolution("resolved" if exact == 1 else "ambiguous", f"Found {exact} exact missense mapping(s).",
                                       provenance, variant, assessments)
        if not assessments:
            status, detail = "outside_selected_reference_index", "No overlapping coding feature is present in the selected transcript index."
        elif all(row.status in CONSEQUENCES for row in assessments):
            status = "unsupported_consequence"
            detail = ", ".join(sorted({row.status for row in assessments}))
        else:
            status = "unresolved_reference_evidence"
            detail = ", ".join(sorted({row.status for row in assessments}))
        return ReferenceResolution(status, detail, provenance, variant, assessments)

    def mappings(self, variant: NormalizedVariant) -> tuple[MissenseMapping, ...]:
        return self.resolve(variant).mappings

    def transcript_context(self, variant: NormalizedVariant) -> tuple[TranscriptContext, ...]:
        """Read pinned local context without repeating consequence resolution.

        The selected index may contain only CDS features; retain that limitation
        instead of treating translation boundaries as measured exon boundaries.
        """
        if normalize(VariantRequest(variant.chromosome,variant.position_1based,variant.ref,variant.alt,
                                    assembly=variant.assembly,organism=variant.organism)) != variant:
            raise ValueError('Normalized reference identity required')
        with self._lock:
            if self._closed: raise RuntimeError('Reference resolver is closed')
            genome=self._chromosome(variant.chromosome)
            if genome.read(variant.position_0based,variant.position_1based)!=variant.ref:
                raise ValueError('Genomic REF differs from context request')
            contexts=[]
            with closing(self._connect()) as db:
                for tid in self._candidate_transcripts(variant.chromosome,variant.position_1based):
                    for pid in self._protein_ids.get(tid,()):
                        audit=self._audit.get((tid,pid))
                        if not audit or audit['status']!=EXACT: continue
                        row=db.execute('SELECT sequence FROM protein WHERE protein_id=? AND transcript_id=?',(pid,tid)).fetchone()
                        if row is None: continue
                        intervals=db.execute('SELECT feature,start_1based,end_1based,strand FROM feature WHERE transcript_id=? AND chromosome=? AND feature IN (\'CDS\',\'exon\') ORDER BY start_1based,end_1based',(tid,variant.chromosome)).fetchall()
                        if any(s!=audit['strand'] for _,_,_,s in intervals):
                            raise ValueError('Reference interval strand differs')
                        cds=tuple(sorted({(a-1,b) for kind,a,b,_ in intervals if kind=='CDS'}))
                        # Source qualification currently proves CDS, not complete
                        # exon annotation. Do not promote arbitrary extra rows.
                        refs=dict(self.provenance,protein_sha256=sha256(row[0].encode()).hexdigest(),
                                  chromosome_fasta_sha256=self._sources[variant.chromosome]['sha256'])
                        contexts.append(TranscriptContext(tid,pid,audit['strand'],(),cds,refs,
                                                         'cds_only' if cds else 'unavailable',row[0]))
            return tuple(contexts)

    @contextmanager
    def _query_sources(self, variant: NormalizedVariant, mapping: MissenseMapping):
        """Read unchanged source bytes for a mapping already resolved by this provider."""
        with self._lock:
            if self._closed:
                raise RuntimeError("Reference resolver is closed")
            genome = self._chromosome(variant.chromosome)
            with closing(self._connect()) as db:
                row = db.execute("SELECT sequence FROM protein WHERE protein_id=? AND transcript_id=?",
                                 (mapping.protein_id, mapping.transcript_id)).fetchone()
            if row is None:
                raise ValueError("Resolved protein is absent from the pinned index")
            yield genome, row[0]

    def resolve_identity(self, variant: NormalizedVariant) -> IdentityResolution:
        result = self.resolve(variant)
        if result.mappings:
            return IdentityResolution(result.mappings, "reference_verified_missense", provenance=result.provenance)
        status = result.status
        if status in {"unsupported", "unsupported_consequence"}:
            return IdentityResolution((), "outside_scoring_scope", "unsupported", unresolved_questions=(result.detail,), provenance=result.provenance)
        if status in {"invalid_request", "reference_mismatch"}:
            return IdentityResolution((), "outside_scoring_scope", "invalid_request" if status == "invalid_request" else "mapping_conflict",
                                      conflicts=(result.detail,), provenance=result.provenance)
        return IdentityResolution((), "reference_unresolved", missing=(status,), unresolved_questions=(result.detail,),
                                  provenance=result.provenance)

    def close(self):
        with self._lock:
            for reader, _ in self._chromosomes.values():
                reader.close()
            self._chromosomes.clear()
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
