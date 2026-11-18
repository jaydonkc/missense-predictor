# Missense predictor

A research prototype for scoring human missense variants and reporting supporting protein and RNA evidence. The pipeline separates the retained DNA and protein score from evidence that is available for interpretation. EVEE serves as a benchmark.

## Run the demo

From the repository root with Python 3.12:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python demo.py
```

Open `demo.ipynb` in Jupyter for the same walkthrough. The example uses synthetic inputs and weights and shows the baseline score, final score, route dispositions, and the evidence used in scoring.

## Research scope

The lab's retained checkpoint has been replayed on 472 archived development alleles using their saved transformed features. The demo uses a synthetic bundle because this preview ships without lab checkpoints, clinical or assay targets, reference databases, or third-party predictor assets. Real scoring requires a separately qualified checkpoint and exact input features. The retained sigmoid score is uncalibrated.

The earlier ASPA/G6PD abundance extension and G6PD growth diagnostic failed their advancement rules. Routed annotations and named effect outputs can be reported separately from scoring. See the CLI options with `python -m evee_next.biology.routed --help`.

## Run with lab artifacts

A lab configuration can reference a qualified cache bundle or checkpoint. Run:

```sh
python -m evee_next.biology.routed request.json --config pipeline.json --policy cache-only
```

The demo script shows the request and configuration shapes. Real artifacts require provenance and version checks.
