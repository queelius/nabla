
# dualr <img src="man/figures/logo.png" align="right" height="139" alt="" />

> Forward-Mode Automatic Differentiation via Dual Numbers

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/dualr)](https://CRAN.R-project.org/package=dualr)
[![R-CMD-check](https://github.com/queelius/dualr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/queelius/dualr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Why automatic differentiation?

Computing derivatives is fundamental to optimization, statistics, and scientific
computing. The three main approaches are:

| | **Finite Differences** | **Symbolic Diff** | **AD (dualr)** |
|---|---|---|---|
| **Accuracy** | O(h) or O(h²) truncation error | Exact | Exact (machine precision) |
| **Cost** | 2p function evaluations | Expression swell for complex f | p forward passes |
| **Control flow** | Works | Breaks on if/for/while | Works through any code |
| **Implementation** | Easy | Requires CAS | Operator overloading |

Finite differences are inaccurate for ill-conditioned problems. Symbolic
differentiation suffers from expression swell and cannot handle control flow.
**Automatic differentiation** gives exact derivatives at machine precision
through any R code — loops, branches, and all.

## How it works

A dual number extends the reals with an infinitesimal ε where ε² = 0:

$$f(x + \varepsilon) = f(x) + f'(x)\,\varepsilon$$

By propagating ε through every arithmetic operation and math function,
`dualr` computes f(x) and f'(x) simultaneously in a single forward pass.
No symbolic manipulation, no finite step sizes — just exact derivatives
via the algebra of dual numbers.

## Installation

```r
# Install from CRAN
install.packages("dualr")

# Or install development version from GitHub
remotes::install_github("queelius/dualr")
```

## Quick start

**Single variable derivative:**

```r
library(dualr)

# Create a dual variable at x = 2 (seeds derivative tracking)
x <- dual_variable(2)

# Evaluate any expression — derivatives propagate automatically
result <- x^3 + sin(x)
value(result)   # f(2) = 8.909
deriv(result)   # f'(2) = 3*4 + cos(2) = 11.584
```

**Multi-parameter score and Hessian (for MLE):**

```r
# Poisson log-likelihood using sufficient statistics
data <- rpois(100, lambda = 3.5)
n <- length(data)
sum_x <- sum(data)

ll <- function(theta) {
  lambda <- theta[1]
  sum_x * log(lambda) - n * lambda
}

# Gradient (score vector) — one call, exact result
score(ll, 3.0)

# Hessian matrix — exact second derivatives via nested duals
hessian(ll, 3.0)

# Observed information = -Hessian
observed_information(ll, 3.0)
```

## Use cases

- **Optimization** — supply exact gradients to `optim()` and `nlminb()`
- **Maximum likelihood estimation** — score vectors, Hessians, and standard errors
- **Sensitivity analysis** — how outputs change with respect to inputs
- **Taylor approximation** — compute polynomial approximations with exact coefficients
- **Curvature analysis** — second-order geometric properties of curves

## Vignettes

<!-- Links work on the pkgdown site; use browseVignettes("dualr") locally -->
- [Introduction to dualr](articles/introduction.html) — dual numbers, arithmetic, composition, comparison with finite differences
- [MLE Workflow](articles/mle-workflow.html) — score, Hessian, observed information, Newton-Raphson on statistical models
- [Higher-Order Derivatives](articles/higher-order.html) — nested duals, curvature, Taylor expansion
- [Optimizer Integration](articles/optimizer-integration.html) — using `score()` and `hessian()` with `optim()` and `nlminb()`

## License

MIT
