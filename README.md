# Human missense predictor

This research release scores human missense variants from prepared molecular features and reports testable assay predictions on supported proteins. EVEE is the clinical benchmark. V1 scores the clinical endpoint. V2 tests a biological hypothesis.

- **V1 clinical score.** The frozen paper-wide heads score a prepared GRCh38 allele. The report records the route and operating cutoff. On 161,579 variants outside historical training and development the model missed 2,689 pathogenic variants versus EVEE's 2,911, and flagged 22,372 benign variants versus EVEE's 12,752. [Architecture, decisions and benchmarks](REPORT.md) · [Full replay receipt](predictor_v1/verification.json)
- **V2 experimental assay profile.** A frozen G6PD model estimates catalytic turnover relative to wild type and normalized human-cell reporter abundance. A lookup covers 9,676 substitutions on one human G6PD reference. The clinical predictor does not consume these experimental estimates. On nine substitutions at held-out residues, V2 correctly predicted both measurements for four, versus three for an abundance-only control. [V2 results and reproducible evaluation](mechanism_v2/README.md)

## Run the predictors

From the repository root with Python 3.12:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m predictor_v1.score predictor_v1/examples/aligned_protein.json
python -m mechanism_v2.report T279A
python -m evee_next.biology.assay_report ENSP00000377192.3 279 T A --compare-alt P --format markdown
```

The first command after installation runs a real clinical checkpoint. New clinical inputs need a qualified reference window and the route-specific molecular features. The original full-context route also requires upstream features and a frozen parent score; a genome coordinate alone cannot reconstruct those inputs. Assay reports identify the endpoint, units and applicable protein reference. The numerical outputs are uncalibrated for clinical use.

The original synthetic walkthrough remains available through `python demo.py` or `demo.ipynb`. A configured lab environment can use `python -m evee_next.biology.routed request.json --config pipeline.json --policy cache-only`.

See the [v1 input contract and worked examples](predictor_v1/README.md) before scoring new variants.
