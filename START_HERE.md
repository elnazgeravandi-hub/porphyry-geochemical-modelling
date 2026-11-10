# Start here

This is a clean public implementation of the stated computational workflow. It contains no original study inputs, private coordinates, personal paths, or hardcoded manuscript outputs. It is not a recovered copy of the original study execution.

## 1. Install

Use Python 3.12. Open a terminal in the extracted project directory:

```bash
python -m venv .venv
```

Activate it on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Or on Linux/macOS:

```bash
source .venv/bin/activate
```

Install the project and all optional dependencies:

```bash
python -m pip install -e ".[report,spatial,dbcv,test]"
```

## 2. Run the synthetic example

```bash
python examples/synthetic_demo/run.py --output results/demo_01 --artificial-validation
```

Use a new output directory for each run. The example generates artificial data locally and runs all nine model branches, the single PC1–C–N fractal reference, and spatial evaluation. It preserves the published optimizer evaluation budgets and can take several minutes or longer.

## 3. Inspect the outputs

- `results/demo_01/verification.json`: successful completion and scope.
- `results/demo_01/model_runs/model_comparison/validation_summary.csv`: model comparison.
- `results/demo_01/model_runs/model_comparison/fairness_audit.json`: paired optimizer conditions.
- `results/demo_01/interpretation/pca/`: rotated loadings, Bartlett scores and metadata.
- `results/demo_01/interpretation/cn/`: shifted scores, fitted thresholds, classes and plot.
- `results/demo_01/exploration/`: capture results and illustrative maps.

These are synthetic outputs, not the manuscript results. Map supports are abstract squares, not real catchments.

## 4. Verify software

```bash
python -m pytest tests -q
```

See `VERIFICATION.md` for the completed verification of this source package.

## 5. Private study execution

Keep authorized inputs outside the public repository, or under the ignored `private_data/` directory. Follow `docs/input_schema.md` and `docs/execution.md`. The study has one contributing catchment per retained sample; fitting remains sample-level. Catchment overlaps are unioned for area evaluation.

Before calling these outputs a numerical reproduction, run the code on the authorized study inputs and compare loadings, scores, thresholds, selected configurations and capture results. If outputs differ, investigate and document the difference; do not insert manuscript values as computed outputs.

## Implementation contract

The exact choices used in this code are defined in `METHOD_SPECIFICATION.md`. They must be described consistently in any manuscript that claims this implementation produced its results. Historical results cannot be verified from a flowchart or figure alone.
