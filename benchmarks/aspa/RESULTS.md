# Can a structure contact distinguish enzyme activity loss from reduced protein amount?

In human aspartoacylase (ASPA), a simple contact rule correctly recognized two of eight substitutions that lost substrate hydrolysis despite preserved protein amount in one study, and one of thirteen in a second study. Both studies measured the chemical conversion of N-acetylaspartate in cell lysates. Neither established turnover per purified active enzyme. The measurements evaluate the proposed protein-mechanism report. They do not evaluate clinical pathogenicity or EVEE.

| Observed in the experiment | Scored substitutions | Correct with abundance alone | Correct with abundance and structure |
|---|---:|---:|---:|
| COS-7, low hydrolysis with preserved amount | 8 | 0 | 2 |
| COS-7, low hydrolysis with reduced amount | 10 | 10 | 10 |
| COS-7, WT-like hydrolysis and amount | 3 | 3 | 3 |
| **COS-7 total** | **21** | **13 (62%)** | **15 (71%)** |
| HEK293, low hydrolysis with preserved amount | 13 | 0 | 1 (8%) |

The structure rule recovered E24D and E285A in COS-7. It predicted reduced amount for five of the other six selective-loss cases. For example, H21G and H21P both lost hydrolysis, but only H21P had reduced protein on the blot. The model predicted both would have reduced amount. At a different residue, it correctly distinguished WT-like C61A and C61S from dual-loss C61W. These same-residue examples illustrate where the rule succeeded and failed within one experiment. [Hershfield et al., Figures 4, 6 and 7](https://pubmed.ncbi.nlm.nih.gov/17391648/)

E285A was the only selective-loss call recovered in HEK293. Eight were called dual loss and four were called WT-like for both endpoints. All thirteen qualifying HEK293 variants had preserved protein on the reported immunoblot, so that study provides no estimate of dual-loss sensitivity. A305E had low hydrolysis in both laboratories, but showed low protein in COS-7 and preserved protein in the HEK293 tagged-construct experiment. The studies differ in cell type, construct and readout. These data cannot identify which difference explains the disagreement. [Hershfield et al.](https://pubmed.ncbi.nlm.nih.gov/17391648/) · [Mendes et al., Figure 3](https://doi.org/10.1002/humu.23181)

## Method and limits

An existing molecular abundance score trained on human G6PD supplied a within-panel percentile for each ASPA substitution. Percentiles below 0.5 predicted reduced amount and hydrolysis. For the others, a WT residue within 5 Å of the catalytic zinc or bound reaction-intermediate analogue predicted selective loss. The structure feature measures WT proximity. It has no substitution-specific mutant coordinates. The cutoff for amount lacks calibration to a WT immunoblot. This hand-specified rule is a diagnostic transfer experiment distinct from the fitted G6PD activity head. [Human ASPA structure 2O4H](https://www.rcsb.org/structure/2O4H)

All 5,843 available molecular predictions were frozen before the paired outcome table was transcribed. Source discovery and an earlier abundance-only study had already exposed some results, so this is retrospective. Five of 26 COS-7 and one of 14 HEK293 records lacked a qualifying endpoint pair under the predeclared criteria. The analysis accounts for them in the exclusions. The observations are few, come from one protein, include overlapping residues and were selected by the experimental authors. The contact rule remains exploratory.

## Reproduce

Run `python3 evaluate.py` beside `predictions.tsv.gz`, `PROTOCOL.md`, `score-lock.json` and `measurements.tsv`. Its output recreates `results.json`, including exclusions, per-variant calls, and both confusion matrices.
