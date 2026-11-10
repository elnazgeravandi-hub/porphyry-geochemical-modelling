# PC1–C–N fractal stage

Consumes saved PCA scores through `--pca-results`; never fits elemental Cu concentrations. Uses shifted PC1 scores, continuous piecewise log–log fitting and data-driven thresholds. Exports dimensionless thresholds, shifted scores, class IDs and a fractal plot. The positive offset is configurable using `--positive-offset` and recorded in metadata. See `docs/implementation_decisions.md` for explicit choices requiring original-study confirmation.

Run `python -m validation_interpretation.cn.run --help`.
