# nabla 0.7.2

* Fixed zero powers at zero and higher derivatives of polynomials at the
  origin, which previously produced `NaN` after the polynomial degree.
* Preserved higher-order exponent derivatives when their first derivative
  vanishes, by checking all nested dual components before simplifying.
* Fixed `cummax()` and `cummin()` on scalar dual numbers to propagate the
  identity derivative. Unsupported vector-valued dual inputs now error.
* Clarified supported operations, floating-point limitations, and the
  dimension/order cost of derivative tensors in the README and metadata.

# nabla 0.7.0

* Added vignette "Higher-Order MLE Analysis" demonstrating third-order
  derivative tensors (`D(f, x, order = 3)`) for computing asymptotic
  skewness of Gamma MLEs, with Monte Carlo validation.
* Renamed package from `dualr` to `nabla`. The S4 class `dualr` retains its
  name (it describes the object type — a dual number in R).
* Added composable total derivative operator `D(f)`:
  - `D(f)` returns the derivative of `f` as a new function.
  - `D(f, x)` evaluates the derivative at `x`.
  - `D(f, x, order = k)` applies `D` k times for k-th order derivative tensors.
  - `D(D(f))` composes naturally for higher-order derivatives.
  - Output tensor shape: each application of `D` appends one n-dimension.
    For `f: R^n -> R`: gradient `(n)`, Hessian `(n,n)`, etc.
    For `f: R^n -> R^m`: Jacobian `(m,n)`, `(m,n,n)`, etc.
* Unified `gradient()`, `hessian()`, and `jacobian()` as thin wrappers
  around `D`, replacing separate seeding strategies with a single composable
  mechanism. This simplifies the codebase at the cost of `O(p)` gradient
  (was `O(1)` passes) and `O(p^2)` Hessian (was `O(p)` passes).
* Removed deprecated second-order functions: `dual2_variable()`,
  `dual2_constant()`, `value2()`, `first_deriv()`, `second_deriv()`,
  `differentiate2()`. Use `dual_variable_n()`, `dual_constant_n()`,
  `deriv_n()`, and `differentiate_n()` instead.
* Removed internal helpers `.make_grad_vector()` and `.make_grad2_vector()`.
* Updated higher-order vignette to use current API and demonstrate `D` operator.

# nabla 0.6.0

* Replaced MLE-specific API with general-purpose derivative functions:
  - `score()` -> `gradient()` — computes the gradient of any scalar-valued
    function (still single-pass via vector-valued derivatives).
  - `hessian()` — unchanged (already mathematically general).
  - `observed_information()` — removed (trivial: just `-hessian()`).
  - `score_and_hessian()` -> `jacobian()` — generalized to compute the
    full m x p Jacobian matrix of any `f: R^p -> R^m`. Accepts functions
    returning lists, numeric vectors, or scalar dualr objects.
* Renamed `R/mle-helpers.R` to `R/derivatives.R`.
* Updated vignettes to use general terminology (gradient/Hessian/Jacobian).

# nabla 0.5.0

* Generalized to arbitrary-order exact derivatives via recursive nesting.
  New API: `dual_variable_n()`, `dual_constant_n()`, `deriv_n()`,
  `differentiate_n()`.
* Removed Rcpp/C++ fast paths — package is now pure R with no compiled code.
  This simplifies installation and aligns with the package's focus on exact
  arbitrary-order derivatives rather than scalar speed.
* Deprecated second-order-only functions (`dual2_variable()`, `dual2_constant()`,
  `value2()`, `first_deriv()`, `second_deriv()`, `differentiate2()`) as thin
  wrappers around the new generalized API.
* Removed benchmarks (speed is not this package's value proposition).

# nabla 0.4.0

* `score()` now computes the full gradient in 1 forward pass (was p passes)
  using vector-valued derivatives, exploiting the `ANY` slots of the `dualr` class.
* `hessian()` now computes the full Hessian in p forward passes (was p(p+1)/2)
  using vector-gradient inner duals with nested outer duals.
* Internal `.is_scalar_dual()` now also checks `length() == 1L` to correctly
  distinguish scalar duals (C++ fast path) from vector-gradient duals (R path).

# nabla 0.3.0

* Added Rcpp-based C++ fast paths for first-order dual arithmetic (`+`, `-`, `*`, `/`, `^`), math (`exp`, `sqrt`, `log`), and `sum`. Provides 3-10x speedup on scalar dual operations while preserving full R fallback for nested (second-order) duals.
* New internal `.is_scalar_dual()` predicate gates C++ vs R paths using `is.double()` on slot contents.
* Added `Rcpp` to `Imports` and `LinkingTo`; package now requires C++ compilation.

# nabla 0.2.0

* Renamed S4 class from `dual` to `dualr` to avoid conflict with base R's `dual` usage.
* Added dedicated `setMethod` dispatches for hot-path arithmetic (`+`, `-`, `*`, `/`, `^`) and math (`exp`, `sqrt`) operations, bypassing group generic overhead.
* Extracted `.dual_min()` / `.dual_max()` internal helpers, deduplicating 6 inline lambdas across `dual-arithmetic.R` and `dual-math.R`.
* Removed dead `switch` branches (`sqrt`, `exp`, `log`) from `Math` group generic that were shadowed by dedicated methods.
* Standardized `sum()` in `Summary` group generic to use `.as_dual()` promotion, consistent with `prod`, `min`, `max`, and `range`.
* Fixed stale `\code{compositional.mle}` reference in `score()` documentation.

# nabla 0.1.0

* Initial CRAN release.
* S4 dual number class with full arithmetic and math function support.
* Nested duals for exact second-order derivatives (`dual2_variable`, `differentiate2`).
* MLE workflow helpers: `score`, `hessian`, `observed_information`, `score_and_hessian`.
* Special functions: `erf`, `erfc`, `beta`, `lbeta`, `psigamma`.
* Four vignettes: introduction, MLE workflow, higher-order derivatives, optimizer integration.
