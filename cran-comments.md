## Update to 0.7.2

This update fixes derivatives of zero powers at the origin, preserves nested
exponent curvature when a first derivative vanishes, and corrects scalar
cumulative extrema derivatives. Documentation now states the supported
operations and numerical/computational limitations explicitly.

## Test environment

* Ubuntu 24.04.3 LTS, R 4.3.3, x86_64-pc-linux-gnu
* R CMD build with vignettes, then R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 2 notes

* Current time could not be verified by the checking environment.
* HTML validation and math rendering checks were skipped because the local
  environment lacks tidy and the R package V8.

All 770 test assertions passed with no failures, warnings, or skips.
Vignette rebuilding, examples, and the PDF manual passed. Current test coverage
is 92.63%. Cross-platform checks have not been rerun for this update.
