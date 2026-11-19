# Standalone human missense predictor

On a matched sample of 10,000 human missense variants, our full-context predictor made 707 errors and EVEE made 786 in the prespecified comparison. The difference was 79 fewer errors, or 0.79 percentage points higher accuracy. In a paired bootstrap that resampled whole genes, the 95% interval was 2.3 to 13.6 fewer errors per 1,000 variants. The study used public retrospective labels.

## Architecture and design

![V1 routes qualified features through either the original full-context model or a paper-wide classifier. The clinical score feeds a variant report and a retrospective comparison with EVEE.](figures/v1-architecture.png)

1. **Resolve the allele.** Align the genomic change to the reference sequence, transcript, codon and amino-acid substitution. An alternative isoform can have a different consequence. Reports identify the transcript and whether it matches the paper.
2. **Score the clinical endpoint.** We compare Evo 2 outputs for the reference and alternate DNA sequences. For the protein, we run frozen Biohub ESMC-600M on both sequences and mask the changed residue to estimate its amino-acid preferences. A second classifier uses up to 104 features, including substitution chemistry, evolutionary constraints and available structure and splice evidence. The model combines its output with the DNA–protein score. The paper-wide route uses three frozen gradient-boosted classifiers per input type. The protein-inclusive heads combine DNA sequence and conservation with substitution chemistry and vertebrate protein alignments. A parallel DNA-only ensemble covers missing protein mappings. EVEE is a comparator. No clinical lookup or assay result enters a new variant's score.
3. **Report molecular effects separately.** For supported proteins, the pipeline reports predictions for specific readouts, including reporter abundance and G6PD-dependent yeast growth. Each output carries its unit and scope. An exploratory alternative-transcript score remains separate from the clinical call.

**Explainability and transparency.** For each call, the workflow records which evidence was available and which fold scores were used. A reader can check whether the residue had structural coverage, inspect the component scores, and reproduce the cutoff decision. Showing that structure or splicing caused disease requires independent biological evidence.

## Matched full-context benchmark against EVEE

| How the cutoff was chosen | Model | Cutoff | Pathogenic missed | Benign flagged | Total errors | Accuracy |
|---|---|---:|---:|---:|---:|---:|
| Fewest errors on the same development set | Our predictor | 0.727886 | 293 | 414 | 707 | 92.93% |
| Fewest errors on the same development set | EVEE | 0.3484 | 310 | 504 | 814 | 91.86% |
| Same numerical cutoff as the paper | Our predictor | 0.5 | 163 | 760 | 923 | 90.77% |
| Same numerical cutoff as the paper | EVEE | 0.5 | 499 | 307 | 806 | 91.94% |

The **prespecified comparison** used our cutoff of 0.7278858715285142 and EVEE’s cutoff of 0.3738. EVEE missed 332 pathogenic variants and flagged 454 benign variants. Our predictor corrected 314 EVEE errors and introduced 235. The cutoffs had been selected on different development sets.

The **later matched-development check** used the same 25,348 development variants to choose each model’s cutoff, counting each error once. Our predictor missed 17 fewer pathogenic variants and flagged 90 fewer benign variants on the 10,000-variant sample. It corrected 340 errors and introduced 233. The gain was 1.07 percentage points in accuracy (95% gene-bootstrap interval 0.46–1.69 points). At 0.5, EVEE made 117 fewer total errors and our predictor missed 336 fewer pathogenic variants. A numerical cutoff of 0.5 gives different trade-offs on the two score scales.

All 10,000 had protein sequence features, 7,400 had a qualified reference structure, and 3,943 had usable burial evidence. The sample was drawn from molecularly verifiable alleles and can include relatives of training proteins. In a separate eligibility-restricted set of 1,956 variants, the prespecified comparison gave 186 versus 206 errors and the paired interval included a tie. [Public scores, development counts, and replay instructions](benchmarks/full_context/README.md)

## Earlier paper-wide benchmark against EVEE

The paper roster has 198,876 variants. The primary comparison used 44,644 pathogenic or likely pathogenic and 116,935 benign or likely benign alleles outside historical TRAIN/DEV. The standalone cutoff (0.43431836, rounded here) was selected on 8,048 historical DEV variants to match EVEE's sensitivity at its project cutoff of 0.1922. Historical investigations had already examined much of this archive.

| Outside-development measure | Standalone | EVEE |
|---|---:|---:|
| Pathogenic variants missed | **2,689** | 2,911 |
| Sensitivity | **93.98%** | 93.48% |
| Benign variants flagged | 22,372 | **12,752** |
| Specificity | 80.87% | **89.09%** |
| AUROC | 0.9465 | **0.9708** |

The sensitivity difference was +0.50 percentage points (95% interval +0.02 to +0.95, resampling genes with the thresholds fixed). The standalone predictor rescued 1,776 pathogenic variants missed by EVEE and missed 1,554 that EVEE caught. [Aggregate counts](benchmarks/missense.json) · [EVEE archive](https://doi.org/10.5281/zenodo.19701997)

**Protein evidence helped on matched alleles.** On 146,123 outside-development alleles with qualified alignments, our DNA-only head missed 4,935 pathogenic variants. Adding protein evidence reduced misses to 2,331 and benign flags from 27,223 to 20,691. EVEE missed 2,775 and flagged 11,962 on the same alleles. Among the 9,875 newly recovered protein mappings, 9,842 involved another transcript. Replacing the original calls with those scores would add 54 net pathogenic misses across the fallback cohort, so the clinical default retains its original calls.

**Same-residue check.** Among 251 retrospective pairs with opposing clinical labels at the same residue, sequence alone ranked 211 correctly, sequence plus the tested reference-structure correction ranked 210, and EVEE ranked 207. The strictly fresh roster yielded zero eligible pairs. The structural correction therefore remains exploratory.
