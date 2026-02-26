---
title: 'nabla: Composable Arbitrary-Order Automatic Differentiation for R'
tags:
  - R
  - automatic differentiation
  - dual numbers
  - numerical computing
  - statistical inference
authors:
  - name: Alexander Towell
    orcid: 0000-0001-6443-9897
    affiliation: 1
affiliations:
  - name: Southern Illinois University Edwardsville
    index: 1
date: 24 February 2026
bibliography: paper.bib
---

# Summary

`nabla` is an R package that computes exact derivatives of arbitrary R functions
through forward-mode automatic differentiation (AD) using dual numbers
[@Baydin2018]. The package provides a single composable operator `D` that
differentiates any R function to any order---exactly, at machine precision,
through loops, branches, and all control flow. Applying `D` once yields the
gradient; applying it twice yields the Hessian; applying it $k$ times yields the
$k$-th order derivative tensor. The convenience functions `gradient()`,
`hessian()`, and `jacobian()` are thin wrappers around `D` for common use cases.

```r
library(nabla)
f <- function(x) x[1]^2 * exp(x[2])
D(f, c(1, 0))              # gradient
D(f, c(1, 0), order = 2)   # Hessian
D(f, c(1, 0), order = 3)   # 2x2x2 third-order tensor
```

`nabla` is implemented in pure R [@RCoreTeam2026] with no compiled code and
depends only on `methods` and `stats`. It is available on CRAN and on GitHub at
<https://github.com/queelius/nabla>.

# Statement of Need

Statistical inference routinely requires exact derivatives. Gradient-based
optimization of log-likelihoods converges faster with exact gradients, and the
observed information matrix---the negative Hessian of the log-likelihood at
the maximum likelihood estimate (MLE)---is the foundation of asymptotic standard
errors and confidence intervals [@Pawitan2001]. Inaccurate Hessians propagate
directly into inaccurate standard errors. Beyond second order, third-order
derivative tensors quantify the asymptotic skewness of MLEs
[@McCullagh1987], a correction that standard Wald intervals miss entirely.

Researchers fitting custom statistical models in R---parametric survival
models, hierarchical likelihoods, or models with masked failure data---often
lack closed-form gradient and Hessian expressions and must rely on numerical
finite differences or hand-coded derivatives. Finite differences introduce
$O(h^2)$ truncation error for central differences that compounds at higher
orders, while hand-coding is tedious and error-prone for models with special
functions like `lgamma`, `digamma`, or `psigamma`. `nabla` eliminates both
problems: users write their objective function in ordinary R, and `D`
returns exact derivatives at any order with no additional effort.

# State of the Field

Several R packages address the problem of computing derivatives.
`numDeriv` [@Gilbert2016] is the most widely used, providing numerical
approximations via Richardson extrapolation. It is reliable for first-order
gradients and adequate Hessians but cannot extend to higher orders without
compounding numerical error, and its accuracy depends on step-size selection.
The `dual` package [@Sartore2022] implements forward-mode AD with dual numbers,
supporting first-order gradients. However, it does not support higher-order
derivatives or provide a composable differentiation operator.
The `madness` package [@Pav2020] provides forward-mode AD for multivariate and
matrix-valued functions, primarily targeting the delta method for approximate
standard errors. It supports first-order derivatives only and uses S3 classes
with matrix storage.
The `salad` package [@Perdry2024] also implements forward-mode AD with S4
classes and supports vectors and matrices, but is limited to first-order
derivatives.
The `autodiffr` package [@Li2018] wraps Julia's `ForwardDiff.jl` to provide
both forward and reverse mode AD in R, but it requires a Julia installation and
is not available on CRAN.
The `Deriv` package [@Clarkson2025] takes a symbolic approach, computing exact
derivatives of R expressions by term rewriting. It handles simple closed-form
functions well but suffers from expression swell on complex models and cannot
differentiate through general control flow (loops, branches, recursion).

`nabla` fills a gap in this landscape by providing arbitrary-order derivatives
through a single composable operator. No existing R package on CRAN supports
exact second-order and higher derivatives via AD without requiring symbolic
algebra or external runtimes. The composable design of `D`---where `D(D(f))`
naturally yields the Hessian, and `D(D(D(f)))` yields the third-order
tensor---distinguishes `nabla` from packages that treat gradient, Hessian, and
higher-order computation as separate problems.

# Software Design

The package defines two S4 classes. The `dualr` class represents a dual number
$a + b\varepsilon$ where $\varepsilon^2 = 0$. Both the `value` and `deriv` slots
accept `ANY` type, which enables recursive nesting: a dual number whose
components are themselves dual numbers. This nesting is the mechanism for
arbitrary-order derivatives---evaluating $f$ on a $k$-fold nested dual extracts
all derivatives up to order $k$ in a single forward pass. The `dual_vector`
class is a list-based container of `dualr` objects supporting `x[i]` indexing,
so that user functions follow the same `x[1]`, `x[2]` convention as R's
`optim()`.

Hot-path arithmetic operators (`+`, `-`, `*`, `/`, `^`) and math functions
(`exp`, `sqrt`, `log`) have direct S4 method dispatch to avoid group generic
overhead. All remaining operations fall through to `Ops`, `Math`, and `Math2`
group generics with `switch`-based dispatch. Special functions---`lgamma`,
`digamma`, `trigamma`, `psigamma`, `erf`, `erfc`, `beta`, and
`lbeta`---propagate derivatives via the chain rule, using the identity
$\frac{d}{dx}\psi^{(n)}(x) = \psi^{(n+1)}(x)$ for the polygamma hierarchy.

The `D` operator is implemented as iterated application of a single-order
differentiation step. Each application runs $p$ forward passes (one per input
dimension), seeding the $j$-th dual vector entry with derivative 1 and all
others with 0. The output at each level is a tensor whose last dimension
indexes the input dimensions. For a scalar function $f: \mathbb{R}^n \to
\mathbb{R}$, `D(f, x)` returns an $n$-vector, `D(f, x, order = 2)` returns an
$n \times n$ matrix, and `D(f, x, order = 3)` returns an $n \times n \times n$
array. Composability is automatic: when the inner `D` returns dual-valued
outputs, the outer `D` propagates through them via the same S4 dispatch.

# Research Impact Statement

`nabla` was developed as part of a broader ecosystem of R packages for
reliability engineering and masked failure data analysis. It provides the
gradient and Hessian infrastructure used by downstream packages including
`algebraic.mle` and `likelihood.model` for Fisher information computation, and
has been used in research on maximum likelihood estimation for series systems
with masked causes of failure. The MLE skewness vignette demonstrates a
concrete application: computing third-order derivative tensors of Gamma
log-likelihoods to quantify asymptotic skewness of MLEs, validated against
Monte Carlo simulation. The package is available on CRAN and is intended to
serve researchers who need exact derivatives for custom statistical models
without the overhead of symbolic algebra systems or external language
runtimes.

# AI Usage Disclosure

Claude Code (Anthropic) was used to assist with drafting this JOSS paper,
including literature search and formatting. The package source code, vignettes,
and all mathematical content were written by the author. All AI-generated text
was reviewed, edited, and verified by the author.

# Acknowledgements

The author thanks the R Core Team for the R language and the CRAN maintainers
for their infrastructure. The S4 method dispatch system, which enables
transparent operator overloading for dual number arithmetic, is central to
this package's design.

# References
