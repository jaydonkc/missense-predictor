# Standalone human missense predictor

On 161,579 paper variants outside historical training and development, the standalone predictor caught 222 more pathogenic or likely pathogenic variants than EVEE and flagged 9,620 more benign or likely benign variants. EVEE had stronger overall discrimination. These results describe a retrospective comparison of frozen scores.

## Architecture and design

1. **Resolve the allele.** Align the genomic change to the reference sequence, transcript, codon and amino-acid substitution. An alternative isoform can have a different consequence. Reports identify the transcript and whether it matches the paper.
2. **Score the clinical endpoint.** The original full-context route uses DNA and protein sequence representations with molecular context where its inputs are qualified. The paper-wide route uses three frozen gradient-boosted classifiers per input type. The protein-inclusive heads combine DNA sequence and conservation with substitution chemistry and vertebrate protein alignments. A parallel DNA-only ensemble covers missing protein mappings. EVEE is a comparator. No clinical lookup or assay result enters a new variant's score.
3. **Report molecular effects separately.** For supported proteins, the pipeline reports predictions for specific readouts, including reporter abundance and G6PD-dependent yeast growth. Each output carries its unit and scope. An exploratory alternative-transcript score remains separate from the clinical call. The V2 experiment below asks whether sequence and structure can distinguish low direct enzyme activity with preserved abundance from low activity accompanied by low abundance.

## Clinical benchmark against EVEE

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

## Did the predictor distinguish biological effects?

The experimental G6PD checkpoint estimates two quantities from sequence and structure: purified enzyme turnover relative to wild type and protein abundance in a human-cell reporter. Assays supplied training targets and held-out measurements. Scoring a new substitution requires only molecular features.

On nine substitutions at five held-out residues, V2 predicted both outcomes correctly for **4/9**, versus **3/9** for an abundance-only control. It identified **4/7** substitutions with low turnover and preserved abundance and missed both substitutions that retained both readouts. One borderline call, R387C, accounted for the gain. Turnover error was essentially identical (0.402 versus 0.401 times wild type). None of the 19 training and evaluation substitutions lost both turnover and abundance, leaving that distinction unverified in G6PD. [G6PD results and assay sources](mechanism_v2/README.md)

A separate ASPA test used reference-structure contacts and a predicted abundance score. In one study it identified **10/10** substitutions with joint loss and only **2/8** with low activity despite preserved protein. In a second study it identified **1/13** of the latter. Selective loss of activity remains the main weakness in these assays. The mechanism reports expose these hypotheses for testing. They have no weight in the clinical score. [ASPA results and assay sources](benchmarks/aspa/RESULTS.md)
