# Contributing to nabla

Thank you for your interest in contributing to nabla! This document
describes how to report bugs, propose changes, and submit pull requests.

## Reporting Bugs

Open an issue at <https://github.com/queelius/nabla/issues> with:

- A minimal reproducible example
- Expected vs actual output
- `sessionInfo()` output

## Proposing Changes

For feature requests or design discussions, open an issue first so we can
align on the approach before you invest time coding.

## Pull Requests

1. Fork the repo and create a feature branch from `main`.
2. Follow the existing coding conventions (see below).
3. Add tests for new functionality in `tests/testthat/`.
4. Run `devtools::test()` and `R CMD check --as-cran` before submitting.
5. Update documentation if you change any exported functions:
   `roxygen2::roxygenise()`.
6. Keep PRs focused — one logical change per PR.

## Coding Style

- **S4 classes and methods** for all dual number types and operators.
- Arithmetic operators (`+`, `-`, `*`, `/`, `^`) and hot-path math
  (`exp`, `sqrt`, `log`) have direct S4 methods; everything else uses
  group generics (`Ops`, `Math`, `Math2`, `Summary`).
- Internal helpers are prefixed with `.` (e.g., `.dual()`, `.as_dual()`).
- Test files mirror source files: `R/dual-math.R` is tested by
  `tests/testthat/test-math.R`.

## Testing

nabla uses testthat edition 3. The standard verification pattern is:

1. Compute the derivative via AD (dual numbers)
2. Compare against the analytical formula
3. Compare against numerical finite differences (`central_difference()`
   from `tests/testthat/helper-numerical.R`)

```r
devtools::test()                    # run all tests
testthat::test_file("tests/testthat/test-math.R")  # run one file
covr::package_coverage()            # check coverage
```

## Code of Conduct

Please note that nabla is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in
this project you agree to abide by its terms.
