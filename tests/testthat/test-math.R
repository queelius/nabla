# Test all math functions against numerical central differences.
# Each test verifies that dual-mode derivative matches the numerical estimate.

tol_deriv <- 1e-6  # tolerance for derivative comparison

# Helper: check a unary function's derivative at a point
check_deriv <- function(f, x, label = deparse(substitute(f))) {
  d <- f(dual(x, 1))
  num <- central_difference(f, x)
  expect_equal(
    deriv(d), num,
    tolerance = tol_deriv,
    label = sprintf("%s'(%g): AD=%g, numerical=%g", label, x, deriv(d), num)
  )
}

# -- Exponential / Logarithm --------------------------------------------------

test_that("exp derivative", {
  check_deriv(exp, 1)
  check_deriv(exp, -2)
  check_deriv(exp, 0)
})

test_that("expm1 derivative", {
  check_deriv(expm1, 0.01)
  check_deriv(expm1, -0.5)
})

test_that("log derivative", {
  check_deriv(log, 1)
  check_deriv(log, 2.5)
  check_deriv(log, 0.1)
})

test_that("log2 derivative", {
  check_deriv(log2, 1)
  check_deriv(log2, 4)
})

test_that("log10 derivative", {
  check_deriv(log10, 1)
  check_deriv(log10, 100)
})

test_that("log1p derivative", {
  check_deriv(log1p, 0.001)
  check_deriv(log1p, 1)
})

test_that("log with base argument", {
  f <- function(x) log(x, base = 3)
  check_deriv(f, 2, "log_base3")
  check_deriv(f, 9, "log_base3")
})

# -- Square root ---------------------------------------------------------------

test_that("sqrt derivative", {
  check_deriv(sqrt, 4)
  check_deriv(sqrt, 0.25)
})

# -- Trigonometric --------------------------------------------------------------

test_that("sin derivative", {
  check_deriv(sin, 0)
  check_deriv(sin, pi / 4)
  check_deriv(sin, pi)
})

test_that("cos derivative", {
  check_deriv(cos, 0)
  check_deriv(cos, pi / 4)
  check_deriv(cos, pi)
})

test_that("tan derivative", {
  check_deriv(tan, 0)
  check_deriv(tan, pi / 4)
  check_deriv(tan, 0.5)
})

# -- Inverse trigonometric ------------------------------------------------------

test_that("asin derivative", {
  check_deriv(asin, 0)
  check_deriv(asin, 0.5)
  check_deriv(asin, -0.5)
})

test_that("acos derivative", {
  check_deriv(acos, 0)
  check_deriv(acos, 0.5)
  check_deriv(acos, -0.5)
})

test_that("atan derivative", {
  check_deriv(atan, 0)
  check_deriv(atan, 1)
  check_deriv(atan, -3)
})

# -- Hyperbolic ----------------------------------------------------------------

test_that("sinh derivative", {
  check_deriv(sinh, 0)
  check_deriv(sinh, 1)
})

test_that("cosh derivative", {
  check_deriv(cosh, 0)
  check_deriv(cosh, 1)
})

test_that("tanh derivative", {
  check_deriv(tanh, 0)
  check_deriv(tanh, 1)
  check_deriv(tanh, -2)
})

# -- Inverse hyperbolic --------------------------------------------------------

test_that("asinh derivative", {
  check_deriv(asinh, 0)
  check_deriv(asinh, 1)
  check_deriv(asinh, -2)
})

test_that("acosh derivative", {
  check_deriv(acosh, 1.5)
  check_deriv(acosh, 3)
})

test_that("atanh derivative", {
  check_deriv(atanh, 0)
  check_deriv(atanh, 0.5)
  check_deriv(atanh, -0.3)
})

# -- Gamma-related -------------------------------------------------------------

test_that("gamma derivative", {
  check_deriv(gamma, 1)
  check_deriv(gamma, 2.5)
  check_deriv(gamma, 0.5)
})

test_that("lgamma derivative", {
  check_deriv(lgamma, 1)
  check_deriv(lgamma, 2.5)
  check_deriv(lgamma, 5)
})

test_that("digamma derivative", {
  check_deriv(digamma, 1)
  check_deriv(digamma, 2.5)
  check_deriv(digamma, 5)
})

test_that("trigamma derivative", {
  check_deriv(trigamma, 1)
  check_deriv(trigamma, 2.5)
  check_deriv(trigamma, 5)
})

# -- Absolute value / sign ------------------------------------------------------

test_that("abs derivative", {
  check_deriv(abs, 3)
  check_deriv(abs, -3)
})

test_that("sign derivative is zero", {
  d <- sign(dual(3, 1))
  expect_equal(value(d), 1)
  expect_equal(deriv(d), 0)
})

# -- Floor / ceiling / round ---------------------------------------------------

test_that("floor derivative is zero", {
  d <- floor(dual(3.7, 1))
  expect_equal(value(d), 3)
  expect_equal(deriv(d), 0)
})

test_that("ceiling derivative is zero", {
  d <- ceiling(dual(3.2, 1))
  expect_equal(value(d), 4)
  expect_equal(deriv(d), 0)
})

# -- atan2 ---------------------------------------------------------------------

test_that("atan2 derivative", {
  # d/dx atan2(x, 2) at x=1 is 2/(1+4) = 2/5
  y <- dual(1, 1)
  r <- atan2(y, 2)
  expect_equal(value(r), atan2(1, 2))
  expect_equal(deriv(r), 2 / 5, tolerance = 1e-10)
})

test_that("atan2 both dual", {
  # d/dx atan2(x, x) = d/dx atan2(1,1) = 0 (constant pi/4)
  # More precisely: (x*1 - x*1)/(x^2+x^2) = 0
  x <- dual(1, 1)
  r <- atan2(x, x)
  expect_equal(value(r), pi / 4)
  expect_equal(deriv(r), 0, tolerance = 1e-10)
})

# -- Compositions ---------------------------------------------------------------

test_that("log(sin(x) + 2) derivative", {
  f <- function(x) log(sin(x) + 2)
  x <- 1.0
  check_deriv(f, x, "log(sin(x)+2)")
})

test_that("exp(cos(x) * x) derivative", {
  f <- function(x) exp(cos(x) * x)
  check_deriv(f, 0.5, "exp(cos(x)*x)")
})

test_that("sqrt(1 + x^2) derivative", {
  f <- function(x) sqrt(1 + x^2)
  check_deriv(f, 2, "sqrt(1+x^2)")
  check_deriv(f, 0, "sqrt(1+x^2)")
})

test_that("lgamma(exp(x)) derivative", {
  f <- function(x) lgamma(exp(x))
  check_deriv(f, 0.5, "lgamma(exp(x))")
})

# -- max / min ------------------------------------------------------------------

test_that("max selects correct branch", {
  a <- dual(3, 1)
  b <- dual(5, 2)
  r <- max(a, b)
  expect_equal(value(r), 5)
  expect_equal(deriv(r), 2)
})

test_that("min selects correct branch", {
  a <- dual(3, 1)
  b <- dual(5, 2)
  r <- min(a, b)
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 1)
})

# -- max / min with 3+ arguments (bug fix regression tests) ------------------

test_that("max with 3 dual arguments selects the largest", {
  a <- dual(3, 10)
  b <- dual(7, 20)
  c <- dual(5, 30)
  r <- max(a, b, c)
  expect_equal(value(r), 7)
  expect_equal(deriv(r), 20)
})

test_that("max with 3 dual arguments when last is largest", {
  a <- dual(1, 10)
  b <- dual(2, 20)
  c <- dual(9, 30)
  r <- max(a, b, c)
  expect_equal(value(r), 9)
  expect_equal(deriv(r), 30)
})

test_that("min with 3 dual arguments selects the smallest", {
  a <- dual(3, 10)
  b <- dual(7, 20)
  c <- dual(5, 30)
  r <- min(a, b, c)
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 10)
})

test_that("min with 3 dual arguments when last is smallest", {
  a <- dual(5, 10)
  b <- dual(3, 20)
  c <- dual(1, 30)
  r <- min(a, b, c)
  expect_equal(value(r), 1)
  expect_equal(deriv(r), 30)
})

test_that("max with mixed dual and numeric args", {
  a <- dual(3, 1)
  r <- max(a, 7, 5)
  expect_true(is_dual(r))
  expect_equal(value(r), 7)
  expect_equal(deriv(r), 0)
})

test_that("min with mixed dual and numeric args", {
  a <- dual(3, 1)
  r <- min(a, 1, 5)
  expect_true(is_dual(r))
  expect_equal(value(r), 1)
  expect_equal(deriv(r), 0)
})

# -- NaN / Inf propagation ----------------------------------------------------

test_that("NaN propagates through dual arithmetic", {
  x <- dual(NaN, 1)
  r <- x + dual(1, 1)
  expect_true(is.nan(value(r)))
})

test_that("Inf propagates through exp of large dual", {
  x <- dual(1000, 1)
  r <- exp(x)
  expect_true(is.infinite(value(r)))
})

test_that("log of negative dual produces NaN", {
  x <- dual(-1, 1)
  r <- log(x)
  expect_true(is.nan(value(r)))
})

# -- cumprod now errors --------------------------------------------------------

test_that("cumprod on dual throws an error", {
  expect_error(cumprod(dual(3, 1)), "cumprod.*not supported")
})
