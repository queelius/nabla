
# nabla <img src="man/figures/logo.png" align="right" height="139" alt="" />

> Arbitrary-order exact derivatives at machine precision

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/nabla)](https://CRAN.R-project.org/package=nabla)
[![R-CMD-check](https://github.com/queelius/nabla/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/queelius/nabla/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`nabla` provides a single composable operator `D` that differentiates any
R function to any order — exactly, at machine precision, through loops,
branches, and all control flow:

```r
library(nabla)
f <- function(x) x[1]^2 * exp(x[2])

D(f, c(1, 0))                # gradient
D(f, c(1, 0), order = 2)     # Hessian
D(f, c(1, 0), order = 3)     # 2×2×2 third-order tensor
D(f, c(1, 0), order = 4)     # 2×2×2×2 fourth-order tensor
```

Each application of `D` adds one dimension to the output. `D(D(f))` gives
the Hessian, `D(D(D(f)))` gives the third-order tensor, and so on —
no limit on order, no loss of precision, no symbolic algebra.

## Why nabla?

| | **Finite Differences** | **Symbolic Diff** | **AD (nabla)** |
|---|---|---|---|
| **Accuracy** | O(h) or O(h²) truncation error | Exact | **Exact (machine precision)** |
| **Higher-order** | Error compounds rapidly | Expression swell | **Composes cleanly to any order** |
| **Control flow** | Works | Breaks on if/for/while | **Works through any code** |

Finite differences lose precision at higher orders (each order multiplies
the error). Symbolic differentiation suffers from expression swell.
`nabla` composes `D` via nested dual numbers — each order is as precise
as the first.

## Installation

```r
# Install from CRAN
install.packages("nabla")

# Or install development version from GitHub
remotes::install_github("queelius/nabla")
```

## The `D` operator

`D` is the core of `nabla`. It differentiates any function `f` and
returns a new function — which can itself be differentiated:

```r
f <- function(x) x[1]^2 * x[2] + sin(x[2])

Df   <- D(f)          # first derivative (function)
DDf  <- D(Df)         # second derivative (function)
DDDf <- D(DDf)        # third derivative (function)

Df(c(3, 4))           # gradient vector
DDf(c(3, 4))          # Hessian matrix
DDDf(c(3, 4))         # 2×2×2 tensor
```

Equivalently, evaluate directly at a point:

```r
D(f, c(3, 4))              # gradient
D(f, c(3, 4), order = 2)   # Hessian
D(f, c(3, 4), order = 3)   # third-order tensor
```

`gradient()`, `hessian()`, and `jacobian()` are convenience wrappers:

```r
gradient(f, c(3, 4))       # == D(f, c(3, 4))
hessian(f, c(3, 4))        # == D(f, c(3, 4), order = 2)
```

## How it works

A dual number extends the reals with an infinitesimal ε where ε² = 0:

$$f(x + \varepsilon) = f(x) + f'(x)\,\varepsilon$$

For higher orders, `nabla` nests dual numbers: a dual whose components are
themselves duals. Each level of nesting extracts one additional order of
derivative — so `D(D(D(f)))` propagates through triply-nested duals to
produce exact third derivatives. This works through `lgamma`, `psigamma`,
trig functions, and all of R's math — no special cases needed.

## Use cases

- **Optimization** — supply exact gradients to `optim()` and `nlminb()`
- **Maximum likelihood** — Hessians for standard errors, third-order tensors
  for asymptotic skewness of MLEs
- **Sensitivity analysis** — how outputs change with respect to inputs
- **Taylor approximation** — exact coefficients to any order
- **Curvature analysis** — second-order geometric properties

## Vignettes

- [Introduction to nabla](https://queelius.github.io/nabla/articles/introduction.html) — dual numbers, arithmetic, composition, comparison with finite differences
- [MLE Workflow](https://queelius.github.io/nabla/articles/mle-workflow.html) — gradient, Hessian, Newton-Raphson on statistical models
- [Higher-Order Derivatives](https://queelius.github.io/nabla/articles/higher-order.html) — the `D` operator, curvature, Taylor expansion
- [Higher-Order MLE Analysis](https://queelius.github.io/nabla/articles/mle-skewness.html) — third-order derivative tensors and asymptotic skewness of MLEs
- [Optimizer Integration](https://queelius.github.io/nabla/articles/optimizer-integration.html) — using `gradient()` and `hessian()` with `optim()` and `nlminb()`

## License

MIT
