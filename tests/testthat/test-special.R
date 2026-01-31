# Tests for special mathematical functions

tol_deriv <- 1e-6

# -- erf / erfc ----------------------------------------------------------------

test_that("erf numeric values are correct", {
  # erf(0) = 0
  expect_equal(erf(0), 0)
  # erf(Inf) = 1
  expect_equal(erf(Inf), 1)
  # erf is odd: erf(-x) = -erf(x)
  expect_equal(erf(-1), -erf(1))
})

test_that("erf derivative is correct", {
  f <- function(x) erf(x)
  for (x in c(0, 0.5, 1, -1, 2)) {
    d <- erf(dual(x, 1))
    num <- central_difference(f, x)
    expect_equal(deriv(d), num, tolerance = tol_deriv,
                 label = sprintf("erf'(%g)", x))
  }
})

test_that("erfc(x) = 1 - erf(x)", {
  for (x in c(0, 0.5, 1, -1, 2)) {
    expect_equal(erfc(x), 1 - erf(x), tolerance = 1e-12)
  }
})

test_that("erfc derivative is correct", {
  f <- function(x) erfc(x)
  for (x in c(0, 0.5, 1, -1)) {
    d <- erfc(dual(x, 1))
    num <- central_difference(f, x)
    expect_equal(deriv(d), num, tolerance = tol_deriv,
                 label = sprintf("erfc'(%g)", x))
  }
})

# -- lbeta ---------------------------------------------------------------------

test_that("lbeta numeric matches base", {
  expect_equal(lbeta(2, 3), base::lbeta(2, 3))
  expect_equal(lbeta(0.5, 0.5), base::lbeta(0.5, 0.5))
})

test_that("lbeta derivative wrt first argument", {
  f <- function(a) lbeta(a, 3)
  for (a in c(0.5, 1, 2, 5)) {
    d <- lbeta(dual(a, 1), 3)
    num <- central_difference(f, a)
    expect_equal(deriv(d), num, tolerance = tol_deriv,
                 label = sprintf("d/da lbeta(%g, 3)", a))
  }
})

test_that("lbeta derivative wrt second argument", {
  f <- function(b) lbeta(2, b)
  for (b in c(0.5, 1, 3, 5)) {
    d <- lbeta(2, dual(b, 1))
    num <- central_difference(f, b)
    expect_equal(deriv(d), num, tolerance = tol_deriv,
                 label = sprintf("d/db lbeta(2, %g)", b))
  }
})

test_that("lbeta derivative wrt both arguments", {
  # d/dx lbeta(x, x) at x=2
  f <- function(x) lbeta(x, x)
  x <- 2
  d <- lbeta(dual(x, 1), dual(x, 1))
  num <- central_difference(f, x)
  expect_equal(deriv(d), num, tolerance = tol_deriv)
})

# -- beta ----------------------------------------------------------------------

test_that("beta numeric matches base", {
  expect_equal(beta(2, 3), base::beta(2, 3))
})

test_that("beta derivative wrt first argument", {
  f <- function(a) beta(a, 3)
  for (a in c(0.5, 1, 2)) {
    d <- beta(dual(a, 1), 3)
    num <- central_difference(f, a)
    expect_equal(deriv(d), num, tolerance = tol_deriv,
                 label = sprintf("d/da beta(%g, 3)", a))
  }
})

test_that("beta derivative wrt second argument", {
  f <- function(b) beta(2, b)
  for (b in c(0.5, 1, 3)) {
    d <- beta(2, dual(b, 1))
    num <- central_difference(f, b)
    expect_equal(deriv(d), num, tolerance = tol_deriv,
                 label = sprintf("d/db beta(2, %g)", b))
  }
})

# -- psigamma ------------------------------------------------------------------

test_that("psigamma(x, 0) matches digamma", {
  for (x in c(1, 2, 3.5)) {
    expect_equal(psigamma(x, 0), digamma(x))
  }
})

test_that("psigamma(x, 1) matches trigamma", {
  for (x in c(1, 2, 3.5)) {
    expect_equal(psigamma(x, 1), trigamma(x))
  }
})

test_that("psigamma dual derivative", {
  # d/dx psigamma(x, 1) = psigamma(x, 2)
  for (x in c(1, 2, 3.5)) {
    d <- psigamma(dual(x, 1), 1L)
    expect_equal(value(d), trigamma(x))
    expect_equal(deriv(d), base::psigamma(x, deriv = 2L), tolerance = 1e-10)
  }
})

# -- lgamma / digamma / trigamma chain ----------------------------------------

test_that("lgamma -> digamma -> trigamma chain rule", {
  # Verify: d/dx lgamma(x) = digamma(x)
  # Verify: d/dx digamma(x) = trigamma(x)
  x <- 2.5
  d_lg <- lgamma(dual(x, 1))
  expect_equal(deriv(d_lg), digamma(x), tolerance = 1e-10)

  d_dg <- digamma(dual(x, 1))
  expect_equal(deriv(d_dg), trigamma(x), tolerance = 1e-10)
})
