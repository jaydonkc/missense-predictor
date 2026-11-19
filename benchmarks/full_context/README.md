# Full-context missense benchmark

Both predictors scored the same 10,000 human missense variants across 3,813 genes. A fixed hash selected this sample from 146,444 molecularly verifiable variants outside historical training and development. The public paper labels include 2,861 pathogenic or likely pathogenic and 7,139 benign or likely benign variants. The original study had previously evaluated 144 of these variants, and related protein families can appear across partitions.

The [original protocol](PROTOCOL.md) preserves the primary comparator and evaluation boundary as recorded before the readout. The same-development-set comparison was added afterward.

## Reproduce the comparison

From the repository root, after installing the project requirements, run:

```sh
python benchmarks/full_context/replay.py
```

The script checks both input files against their SHA-256 hashes. It recomputes the cutoffs from development counts, reconstructs each predictor score from its eligible folds, compares archived EVEE scores with the same public labels, and checks the results against the recorded counts and paired intervals. It runs on CPU. [Source archive for EVEE and clinical labels](https://doi.org/10.5281/zenodo.19701997)

- `scored_variants.tsv.gz` contains the locked allele ID, gene ID, eligible folds, three frozen fold scores, combined score, archived EVEE score, and public binary endpoint. The predictor used the single eligible held fold for 85 alleles and averaged all three for 9,915. Labels and EVEE scores enter only during evaluation.
- `evidence_availability.tsv.gz` records which evidence was available for each variant, how many context features were observed, and the names of missing fields. Join it to the score table by variant ID. A missing structure or splice observation has no implied benign effect.
- `development_counts.tsv.gz` records the number of positive and negative labels at each distinct score for each predictor on the same 25,348 development alleles. This is sufficient to reproduce both minimum-error cutoffs. A score equal to the cutoff is called pathogenic. The selection rule minimizes false positives plus false negatives and breaks ties by proximity to 0.5, then by the higher cutoff.
- `results.json` includes the original prespecified comparison, the later same-development-set check and the common 0.5 cutoff. Intervals use 5,000 paired whole-gene bootstrap draws with fixed cutoffs and seed 1729. Positive interval values mean fewer errors for our predictor.
- `manifest.json` identifies the source archive, locked score and outcome artifacts, and the public extracts by hash.

## How the DGX run produced the scores

1. Validate the GRCh38 substitution and the residue on its versioned human reference protein. Freeze the allele roster and resolve the eligible training folds before joining the paper labels.
2. Run Evo 2 on reference and alternate DNA windows in both orientations. Extract and project the sequence representations. Run frozen Biohub ESMC-600M on the reference and substituted protein windows. A separate masked pass scores all 20 amino acids at the changed position. Long proteins use windows of at most 1,022 residues.
3. Build the parent score for each fold from the DNA and protein representations and their qualified context. Assemble 104 ordered features for the context classifier, including substitution chemistry, conservation and available reference structure and splice evidence. Source failures are tracked separately from verified gaps in optional evidence.
4. For each fold, combine the frozen parent score with the matching context head. Apply the preassigned fold rule to obtain one score per variant. Freeze the scores before joining public outcomes. The independent lab replay covered 30,000 fold outputs. The run reused 147 cases. For 9,853 newly processed variants, DNA extraction took 47.49 B200 GPU-hours and the protein jobs took 1.56 GPU-hours. CPU annotation and retrieval are additional.

The public benchmark starts with the frozen outputs from that workflow. To repeat the final model calculation with the lab's prepared DGX intermediates, run:

```sh
python benchmarks/full_context/replay.py --prepared /path/to/intermediates.npz
```

The archive must contain the locked variant IDs (`ids`), the 10,000 × 104 feature matrix (`x`, in the model manifest's order), and the 10,000 × 3 parent scores (`adjusted_parent`, in fold order). This option loads the published checkpoint, recomputes all three fold scores and checks them against the public records within an absolute tolerance of 1e-12. Regenerating the molecular inputs from raw variants requires the lab's upstream models, projections, indexed references and caches. The public CPU replay requires none of those assets.
