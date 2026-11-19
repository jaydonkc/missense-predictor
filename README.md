# Human missense predictor

This research release scores human missense variants from prepared molecular features and reports testable assay predictions on supported proteins. EVEE is the clinical benchmark.

**Full-context benchmark.** On the same 10,000 human missense variants outside historical training and development, our predictor made 707 errors versus EVEE’s 786 at the prespecified development-selected cutoffs. A later comparison choosing both cutoffs on the same development variants gave 707 versus 814. At the common paper cutoff of 0.5, EVEE made 806 errors versus our predictor’s 923. [Four-row comparison and architecture](REPORT.md) · [Data, DGX workflow and CPU replay](benchmarks/full_context/README.md)

- **Paper-wide scoring.** A separate, earlier benchmark used fallback routes for all 161,579 outside-development variants. The report retains those results and identifies the inputs each route needs.
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
