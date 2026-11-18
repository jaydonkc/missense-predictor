# Human G6PD turnover and abundance

This research checkpoint predicts two measurements for each human G6PD substitution on ENSP00000377192.3. One estimate is purified enzyme turnover (kcat) relative to matched wild type. The other is steady-state abundance in a human-cell reporter. A low kcat with preserved abundance is consistent with catalytic impairment in the purified preparation. Reduced reporter abundance can reflect several cellular processes. Neither measurement proves a mechanism of disease in red blood cells.

## Run a prediction or reproduce the evaluation

The commands run from this directory using Python 3 and its standard library.

```sh
python3 predict.py --variant R387C
python3 report.py R387C
python3 report.py R387C --json
python3 evaluate.py
```

Each prediction uses frozen molecular scores. The inference commands open only the model, prediction lock, and molecular prediction table. The evaluation command independently joins the released turnover and reporter extracts and verifies the per-substitution and aggregate results. The fitted heads can also score four supplied molecular features using `predict.py --features`. These values must use the original SaProt compatibility, ESMC likelihood, ThermoMPNN stability-change and ASPA-trained abundance percentile scales. The package includes 9,676 scored substitutions with qualified molecular features. Substitutions outside this set need independently qualified inputs. The reference training script documents the fit using NumPy and the historical files identified by SHA-256 in the score lock.

## Held-residue results

The turnover head was fit on 10 substitutions at 10 residues. The frozen ASPA-trained abundance head was mapped onto G6PD reporter units using 190 substitutions at those same residues. Four external biochemical papers supplied 9 evaluable substitutions at 5 other residues. All 9 had preserved reporter abundance. The complete predicted assay profile was correct for **4/9** versus **3/9** using a control that inferred turnover from the molecular abundance score alone. Turnover mean absolute error was **0.402 versus 0.401** times wild-type turnover. Reporter mean absolute error was **0.279** assay units for both.

| Measured profile | Variants | Both endpoints correct | Control correct |
|---|---:|---:|---:|
| Low turnover, preserved abundance | 7 | 4 | 3 |
| Low turnover, low abundance | 0 | No observations | No observations |
| Retained turnover, preserved abundance | 2 | 0 | 0 |
| Retained turnover, low abundance | 0 | No observations | No observations |

The new head correctly called all three C385 substitutions and R387C. It missed all three G410 substitutions on both endpoints. For G410A the measured turnover was 0.162 times wild type and reporter abundance was 1.345. The prediction was 0.524 and 0.731, crossing each decision boundary in the wrong direction. At R387C, observed turnover was 0.113 and the predicted turnover of 0.475 narrowly crossed the 0.50 cutoff. The control predicted 0.511. That one call accounts for the entire difference in exact profiles. The model also missed the complete profiles of both near-WT-turnover variants, R136C and M159I.

Three of seven low-turnover examples share C385 and three share G410. Nine observations therefore give limited information about new sites. The sample contains no measured dual-loss substitutions, so it cannot establish whether the two low-turnover profiles can be distinguished. No clinical decision uses this model. Frozen boundaries, data split and the retrospective nature of this evaluation are recorded in [the protocol](PROTOCOL.md).

## Measurement sources and exclusions

Purified turnover targets were divided by wild-type turnover in the same assay. Training used the public measurements of [Gómez-Manzo et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4691071/), [Boonyuen et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5625996/) and [Boonyuen et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8692357/). Human-cell reporter abundance came from the Geck et al. G6PD abundance dataset. The held extract includes its individual MaveDB score identifiers. [MaveDB 00001266-b-1](https://www.mavedb.org/score-sets/urn:mavedb:00001266-b-1)

The four assessment panels included 25 constructs in total. Four wild-type controls supplied denominators, four were compound substitutions, three overlapped a training residue, two had conflicting test-construct numbering and three lacked a finite kcat. All nine eligible substitutions had measured reporter abundance.

- Martínez-Rosas et al., Table 2, supplied R387C and matched wild type. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC7215812/)
- Boonyuen et al., Table 3, supplied R136C, M159I and matched wild type. Table 3 reports H27R and V298M while Table 4 lists the corresponding sites as H32R and V291M. Both entries were excluded. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC10835850/)
- Gálvez-Ramírez et al., Table 1, supplied C385R, C385G, C385W and matched wild type. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC11641639/)
- Hernández-Ochoa et al., Table 1, supplied G410A, G410C, G410D and matched wild type. Three P409 variants had activity below the assay's quantifiable range and lacked numeric turnover. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC12429737/)

The molecular score lock predates the held reporter join. Some source activity values were visible during the literature search. The panels cover a single protein and reflect the researchers' variant choices. This constitutes retrospective, residue-held assessment. It supplies neither a blinded prospective test nor a clinical comparison to EVEE.

A secondary count checked the ten training substitutions with their original historical turnover measurements and the same frozen cutoffs. Five had low turnover with preserved abundance, four retained both, and one (G163S) retained turnover with low abundance. Across the ten training and nine held substitutions there were **zero dual-loss observations**. The training cases do not contribute to the 4/9 evaluation metric. [Counts and individual training profiles](training-profile-counts.json)
