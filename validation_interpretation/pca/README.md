# PCA stage

Pearson correlation PCA of unscaled clr coordinates, three retained components, Kaiser-normalized Varimax, and Bartlett WLS scoring. Outputs include PC1–PC3 loadings, scores, variance and score coefficients. Only PC1 feeds the subsequent C–N fractal reference. PCA scores alone are not a separately evaluated reference. See `docs/implementation_decisions.md` for the exact formula and study-confirmation requirements.

Run `python -m validation_interpretation.pca.run --help`.
