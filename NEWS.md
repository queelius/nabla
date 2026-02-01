# dualr 0.2.0

* Renamed S4 class from `dual` to `dualr` to avoid conflict with base R's `dual` usage.
* Added dedicated `setMethod` dispatches for hot-path arithmetic (`+`, `-`, `*`, `/`, `^`) and math (`exp`, `sqrt`) operations, bypassing group generic overhead.
* Extracted `.dual_min()` / `.dual_max()` internal helpers, deduplicating 6 inline lambdas across `dual-arithmetic.R` and `dual-math.R`.
* Removed dead `switch` branches (`sqrt`, `exp`, `log`) from `Math` group generic that were shadowed by dedicated methods.
* Standardized `sum()` in `Summary` group generic to use `.as_dual()` promotion, consistent with `prod`, `min`, `max`, and `range`.
* Fixed stale `\code{compositional.mle}` reference in `score()` documentation.

# dualr 0.1.0

* Initial CRAN release.
* S4 dual number class with full arithmetic and math function support.
* Nested duals for exact second-order derivatives (`dual2_variable`, `differentiate2`).
* MLE workflow helpers: `score`, `hessian`, `observed_information`, `score_and_hessian`.
* Special functions: `erf`, `erfc`, `beta`, `lbeta`, `psigamma`.
* Four vignettes: introduction, MLE workflow, higher-order derivatives, optimizer integration.
