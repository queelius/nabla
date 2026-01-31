# Tests for second-order derivatives via nested duals

tol_deriv <- 1e-5  # second-order needs slightly more tolerance

# -- dual2 constructors --------------------------------------------------------

test_that("dual2_variable structure", {
  x <- dual2_variable(3)
  expect_true(is_dual(x))
  expect_true(is_dual(value(x)))
  expect_true(is_dual(deriv(x)))

  expect_equal(value2(x), 3)
  expect_equal(first_deriv(x), 1)
  expect_equal(second_deriv(x), 0)
})

test_that("dual2_constant structure", {
  x <- dual2_constant(5)
  expect_equal(value2(x), 5)
  expect_equal(first_deriv(x), 0)
  expect_equal(second_deriv(x), 0)
})

# -- Second derivatives of basic functions ------------------------------------

test_that("x^2: f''=2", {
  x <- dual2_variable(3)
  r <- x^2
  expect_equal(value2(r), 9)
  expect_equal(first_deriv(r), 6)
  expect_equal(second_deriv(r), 2)
})

test_that("x^3: f''=6x", {
  x <- dual2_variable(2)
  r <- x^3
  expect_equal(value2(r), 8)
  expect_equal(first_deriv(r), 12)  # 3*4 = 12
  expect_equal(second_deriv(r), 12) # 6*2 = 12
})

test_that("exp(x): f''=exp(x)", {
  x_val <- 1.5
  x <- dual2_variable(x_val)
  r <- exp(x)
  expect_equal(value2(r), exp(x_val))
  expect_equal(first_deriv(r), exp(x_val), tolerance = 1e-10)
  expect_equal(second_deriv(r), exp(x_val), tolerance = 1e-10)
})

test_that("sin(x): f''=-sin(x)", {
  x_val <- pi / 3
  x <- dual2_variable(x_val)
  r <- sin(x)
  expect_equal(value2(r), sin(x_val))
  expect_equal(first_deriv(r), cos(x_val), tolerance = 1e-10)
  expect_equal(second_deriv(r), -sin(x_val), tolerance = 1e-10)
})

test_that("cos(x): f''=-cos(x)", {
  x_val <- pi / 4
  x <- dual2_variable(x_val)
  r <- cos(x)
  expect_equal(value2(r), cos(x_val))
  expect_equal(first_deriv(r), -sin(x_val), tolerance = 1e-10)
  expect_equal(second_deriv(r), -cos(x_val), tolerance = 1e-10)
})

test_that("log(x): f''=-1/x^2", {
  x_val <- 2
  x <- dual2_variable(x_val)
  r <- log(x)
  expect_equal(value2(r), log(x_val))
  expect_equal(first_deriv(r), 1/x_val, tolerance = 1e-10)
  expect_equal(second_deriv(r), -1/x_val^2, tolerance = 1e-10)
})

test_that("1/x: f''=2/x^3", {
  x_val <- 3
  x <- dual2_variable(x_val)
  r <- 1 / x
  expect_equal(value2(r), 1/x_val)
  expect_equal(first_deriv(r), -1/x_val^2, tolerance = 1e-10)
  expect_equal(second_deriv(r), 2/x_val^3, tolerance = 1e-10)
})

test_that("sqrt(x): f''=-1/(4*x^(3/2))", {
  x_val <- 4
  x <- dual2_variable(x_val)
  r <- sqrt(x)
  expect_equal(value2(r), 2)
  expect_equal(first_deriv(r), 0.25, tolerance = 1e-10)
  expect_equal(second_deriv(r), -1/(4*x_val^(3/2)), tolerance = 1e-10)
})

# -- Second derivatives against numerical -------------------------------------

test_that("second derivatives match numerical for compositions", {
  fns <- list(
    list(f = function(x) x^4, name = "x^4"),
    list(f = function(x) exp(sin(x)), name = "exp(sin(x))"),
    list(f = function(x) log(1 + x^2), name = "log(1+x^2)"),
    list(f = function(x) x * exp(-x), name = "x*exp(-x)"),
    list(f = function(x) 1 / (1 + exp(-x)), name = "sigmoid")
  )

  for (fn in fns) {
    x_val <- 1.0
    result <- differentiate2(fn$f, x_val)
    num_second <- numerical_second_deriv(fn$f, x_val)
    expect_equal(
      result$second, num_second,
      tolerance = tol_deriv,
      label = sprintf("f''(%g) for %s", x_val, fn$name)
    )
  }
})

# -- differentiate2 helper ------------------------------------------------------

test_that("differentiate2 returns correct structure", {
  result <- differentiate2(sin, pi/4)
  expect_true(is.list(result))
  expect_named(result, c("value", "first", "second"))
  expect_equal(result$value, sin(pi/4))
  expect_equal(result$first, cos(pi/4), tolerance = 1e-10)
  expect_equal(result$second, -sin(pi/4), tolerance = 1e-10)
})

test_that("differentiate2 with lgamma", {
  x_val <- 2.5
  result <- differentiate2(lgamma, x_val)
  expect_equal(result$value, lgamma(x_val))
  expect_equal(result$first, digamma(x_val), tolerance = 1e-10)
  expect_equal(result$second, trigamma(x_val), tolerance = 1e-10)
})

# -- Nested dual power rule (.is_zero fix) ------------------------------------

test_that("x^3 second derivative via nested duals uses power rule path", {
  # f(x) = x^3, f''(x) = 6x
  # At x=2: f''=12
  x <- dual2_variable(2)
  r <- x^3
  expect_equal(second_deriv(r), 12)
})

test_that("2^x second derivative via nested duals", {
  # f(x) = 2^x, f'(x) = 2^x * log(2), f''(x) = 2^x * log(2)^2
  # At x=3: f''= 8 * log(2)^2
  x <- dual2_variable(3)
  r <- 2^x
  expect_equal(second_deriv(r), 8 * log(2)^2, tolerance = 1e-10)
})
