# Frozen v1 predictor

From the repository root, install the dependencies in the root requirements file. The following commands score real molecular features with the frozen checkpoints used in the paper benchmark.

```sh
python -m predictor_v1.score predictor_v1/examples/aligned_protein.json
python -m predictor_v1.score predictor_v1/examples/dna.json
python -m predictor_v1.score predictor_v1/examples/protein_chemistry.json
python -m predictor_v1.score predictor_v1/examples/full_context.json
python -m unittest predictor_v1.test_score
```

The returned score is uncalibrated. A score at or above 0.4343183620275828 is screen-positive at the development-selected threshold. The executable, source checkpoint hashes, weight checksum and expected numerical outputs are included. The four variants in the examples are public GRCh38 benchmark alleles. [Frozen model provenance](assets/manifest.json) · [Example scores](examples/expected.json) · [Full-roster replay](verification.json)

## Input contract

Each input specifies a primary GRCh38 chromosome, zero-based genomic position, and literal reference and alternate DNA bases in forward genomic orientation. The caller supplies molecular features from a verified mapping. The scorer checks lengths, identifiers, values and the central base, but does not fetch an alignment or establish the disease-relevant isoform.

- **DNA only.** Supply the 11-base forward genomic window, centered on the SNV, and the UCSC hg38 phyloP100way score. Use JSON `null` for a missing track or window. The encoder orients the substitution to the pyrimidine strand and includes DNA flanks, GC fraction, CpG status and the substituted bases.
- **Protein chemistry.** Also supply the versioned Ensembl transcript and protein, one-based residue number and two amino acids. Verify the genome, codon, reference peptide, and molecular consequence before scoring. When the alignment is absent, the protein head uses chemistry and its learned missing-feature branches.
- **Aligned protein.** Add the complete twelve-field UCSC hg38 vertebrate protein alignment profile shown in [the example](examples/aligned_protein.json). It summarizes effective species count, weighted reference and replacement frequencies, entropy, reference-to-replacement log odds, BLOSUM-adjusted log odds, nearest callable evolutionary distance, nearest replacement support, clade coverage and primate support. Profiles use the original branch-length weights and at least eight callable nonhuman species. A partial or arbitrary alignment is outside this interface. Both broad heads average three fixed gradient-boosting classifiers for new variants.
- **Full context.** Supply a previously qualified molecular feature vector, original parent score and original held fold as in [the example](examples/full_context.json). The scorer applies the midpoint rule to that fold's model. Assigning folds or substituting another genomic model's score for the frozen parent would change this prediction. Upstream reference databases, embedding extractors and the full-roster input cache are external. The compact release accepts prepared features and cannot derive this route from a bare genomic coordinate.

The single-variant CLI rejects unspecified replacements. The historical benchmark counted seven ambiguous records via a separate predeclared combination of literal-allele scores. Its numerical replay used that historical policy.

## Testable assay report

```sh
python -m evee_next.biology.assay_report ENSP00000377192.3 279 T A --compare-alt P --format markdown
```

Omit `--compare-alt` to report one substitution. The report retrieves frozen model predictions by exact versioned protein, position and substitution. It returns an explicitly missing result when an otherwise supported site lacks a qualified endpoint. Supported proteins for abundance are ASPA, G6PD and SOD1. Yeast growth is covered for human G6PD. The frozen packet contains 18,151 abundance predictions and 9,676 growth predictions. No assay results enter a request. The packet stores predictions computed for these three proteins; the protein-only reporter does not emit a genomic clinical score.

For G6PD T279A and T279P, the model predicts abundance percentiles 75.6 and 70.3 and normalized yeast growth 0.927 and 0.814. A growth effect can reflect abundance, function or both. See the study report for evaluation boundaries and the direct-activity experiment.

The clinical comparison uses archived [EVEE predictions](https://zenodo.org/records/19701997). Neither clinical labels nor EVEE predictions are packaged as model inputs.

The [genomic pair report](examples/genomic-G6PD-report.md) and [alternative-transcript example](examples/alternative-transcript-report.md) show how the research pipeline recorded original scores, missing information and exploratory evidence. The genomic reporting CLI supports a separately supplied reference configuration. Its optional `--paper-study` lookup requires the matching external export and its checksum receipts.
