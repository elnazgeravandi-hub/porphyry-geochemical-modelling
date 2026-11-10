# Public implementation: method specification

This document states executable conventions. It does not assert that unspecified historical study steps used these conventions.

| Stage | Implementation |
|---|---|
| Input | Six positive or explicitly censored concentrations; unique sample IDs; supplied units |
| Quality control | Explicit censoring flags, finite positive measured values, input alignment and units |
| Censoring | Independently implemented ALR/RF upper-truncation procedure; 100 trees, depth 8, seed 42; OOB RMSE scale; synchronous updates; convergence tolerance 1e-5; maximum 100 iterations |
| Transformation | Unit conversion to ppm followed by six-part clr |
| Clustering scaling | Column median and IQR, fitted once and shared within each comparison |
| PCA | Pearson correlation of unscaled clr variables; three retained leading nonzero components; Kaiser-normalized Varimax |
| Scoring | Bartlett weighted least squares with residual variances 1 minus retained communality; reject nonpositive residual variances and ill-conditioned systems |
| Component convention | Keep rotation order; first rotated component supplies PC1; signs oriented to positive Cu |
| Positive scores | s = PC1 - min(PC1) + offset; default offset 1.0; offset and resulting shift are recorded |
| Fractal fit | Distinct positive score thresholds and cumulative sample counts; continuous log10 hinge regressions with one to five segments; minimum 20 distinct thresholds per segment |
| Breakpoints | Numerical differential-evolution search; BIC includes fitted breakpoint positions; sequential additional-segment acceptance requires BIC reduction of at least 6 |
| Classification | Fitted score breakpoints; exact breakpoint ties enter the upper class; no fixed study thresholds |
| Reference | Only the combined PC1–C–N fractal classification; no separately evaluated PCA-score or Cu-fractal reference |
| Clustering | Three baselines and six optimized branches; manuscript domains, constraints and budgets are explicit in configuration and implementation documentation |
| Selection | Feasible unique partitions pooled across seeds; full-data reevaluation; nondominated filtering and common normalized distance-to-ideal selection |
| Spatial units | Sample-linked contributing catchments, clipped and unioned for target area; geometric overlap is not counted twice |
| Occurrence capture | Unique eligible locations in whole-group cumulative targets at common 5%, 10%, and 20% nominal ceilings; actual selected area also reported |
| Study privacy | No original concentrations, identifiers, coordinates, catchments or figures included |
| Demonstration | Independent synthetic generator; abstract spatial supports; execution demonstration only |

## Important numerical scope

Three retained components are a specified setting, not an inferred retention criterion. Bartlett coefficients, score-shifting constants, rotation ordering and BIC breakpoint selection must be checked against the study if exact numerical reproduction is claimed. The current code is transparent about these choices instead of inferring them from a published figure.

Synthetic class counts and thresholds are data dependent. The number of segments is not forced to reproduce the original figure. The optimizer search and breakpoint search are numerical procedures without a claim of a globally optimal solution.

See `docs/implementation_decisions.md` for equations, decoding, numerical details and additional conventions.
