# Tests for base R interoperability: sum(), prod(), c(), is.numeric()
# These verify that duals behave like regular R numbers in common patterns.

# -- sum() ---------------------------------------------------------------------

test_that("sum of two duals", {
  a <- dual(3, 1)
  b <- dual(5, 2)
  r <- sum(a, b)
  expect_true(is_dual(r))
  expect_equal(value(r), 8)
  expect_equal(deriv(r), 3)
})

test_that("sum of dual and numeric", {
  a <- dual(3, 1)
  r <- sum(a, 10)
  expect_true(is_dual(r))
  expect_equal(value(r), 13)
  expect_equal(deriv(r), 1)
})

test_that("sum of three duals", {
  a <- dual(1, 1)
  b <- dual(2, 0)
  c_val <- dual(3, 0)
  r <- sum(a, b, c_val)
  expect_equal(value(r), 6)
  expect_equal(deriv(r), 1)
})

test_that("sum of single dual", {
  a <- dual(5, 3)
  r <- sum(a)
  expect_true(is_dual(r))
  expect_equal(value(r), 5)
  expect_equal(deriv(r), 3)
})

# -- prod() -------------------------------------------------------------------

test_that("prod of two duals (product rule)", {
  a <- dual(3, 1)
  b <- dual(5, 2)
  r <- prod(a, b)
  expect_true(is_dual(r))
  expect_equal(value(r), 15)
  # d/dx [3*5] with da=1, db=2: 1*5 + 3*2 = 11
  expect_equal(deriv(r), 11)
})

test_that("prod of dual and numeric", {
  a <- dual(4, 1)
  r <- prod(a, 3)
  expect_true(is_dual(r))
  expect_equal(value(r), 12)
  expect_equal(deriv(r), 3)
})

test_that("prod of three duals", {
  # d/dx [x * 2 * 4] = 8 at any x (with dx=1)
  x <- dual(3, 1)
  a <- dual(2, 0)
  b <- dual(4, 0)
  r <- prod(x, a, b)
  expect_equal(value(r), 24)
  expect_equal(deriv(r), 8)
})

# -- c() ----------------------------------------------------------------------

test_that("c() creates dual_vector from duals", {
  a <- dual(1, 0)
  b <- dual(2, 1)
  r <- c(a, b)
  expect_true(is(r, "dual_vector"))
  expect_equal(length(r), 2)
  expect_equal(value(r[1]), 1)
  expect_equal(value(r[2]), 2)
  expect_equal(deriv(r[2]), 1)
})

test_that("c() promotes numeric to dual constant", {
  a <- dual(1, 1)
  r <- c(a, 5)
  expect_true(is(r, "dual_vector"))
  expect_equal(length(r), 2)
  expect_equal(value(r[2]), 5)
  expect_equal(deriv(r[2]), 0)
})

test_that("c() with three duals", {
  a <- dual(1, 0)
  b <- dual(2, 1)
  d <- dual(3, 0)
  r <- c(a, b, d)
  expect_equal(length(r), 3)
  expect_equal(value(r[3]), 3)
})

# -- is.numeric() -------------------------------------------------------------

test_that("is.numeric returns TRUE for dual", {
  expect_true(is.numeric(dual(3, 1)))
  expect_true(is.numeric(dual_variable(5)))
  expect_true(is.numeric(dual_constant(0)))
})

# -- if/else with dual comparison (should work) --------------------------------

test_that("if-else with dual comparison works", {
  x <- dual(3, 1)
  result <- if (x > 0) x^2 else -x
  expect_true(is_dual(result))
  expect_equal(value(result), 9)
  expect_equal(deriv(result), 6)
})

# -- Reduce patterns -----------------------------------------------------------

test_that("Reduce('+') works as sum alternative", {
  duals <- list(dual(1, 0.5), dual(2, 1), dual(3, 0.5))
  r <- Reduce("+", duals)
  expect_equal(value(r), 6)
  expect_equal(deriv(r), 2)
})

test_that("Reduce('*') works as prod alternative", {
  duals <- list(dual(2, 1), dual(3, 0), dual(4, 0))
  r <- Reduce("*", duals)
  expect_equal(value(r), 24)
  expect_equal(deriv(r), 12)
})

# -- sapply with accessors -----------------------------------------------------

test_that("sapply(duals, value) extracts values", {
  duals <- list(dual(1, 0), dual(2, 1), dual(3, 0))
  vals <- sapply(duals, value)
  expect_equal(vals, c(1, 2, 3))
})

test_that("sapply(duals, deriv) extracts derivatives", {
  duals <- list(dual(1, 0), dual(2, 1), dual(3, 2))
  derivs <- sapply(duals, deriv)
  expect_equal(derivs, c(0, 1, 2))
})

# -- For loop accumulation -----------------------------------------------------

test_that("for loop accumulation with duals", {
  x <- dual_variable(2)
  accum <- dual(0, 0)
  for (i in 1:3) {
    accum <- accum + x^i
  }
  # x + x^2 + x^3 at x=2: 2+4+8 = 14
  expect_equal(value(accum), 14)
  # 1 + 2x + 3x^2 at x=2: 1+4+12 = 17
  expect_equal(deriv(accum), 17)
})

# -- MLE-style: sum over data with dual parameter ----------------------------

test_that("element-wise residual sum via lapply + Reduce", {
  data_vals <- c(1, 2, 3, 4, 5)
  mu <- dual_variable(3)
  residuals <- lapply(data_vals, function(xi) (xi - mu)^2)
  ss <- Reduce("+", residuals)
  ll <- -0.5 * ss
  # sum((x-3)^2) = 4+1+0+1+4 = 10
  expect_equal(value(ll), -5)
  # d/dmu [-0.5 * sum((xi-mu)^2)] = sum(xi - mu) = -3+(-1)+0+1+2 = 0-1 wait
  # Actually: d/dmu [-0.5 * sum((xi-mu)^2)] = sum(xi - mu) = (1-3)+(2-3)+(3-3)+(4-3)+(5-3)
  # = -2 + -1 + 0 + 1 + 2 = 0
  expect_equal(deriv(ll), 0)
})

# -- Map / vapply patterns -----------------------------------------------------

test_that("Map with two dual lists", {
  xs <- list(dual(1, 1), dual(2, 0))
  ys <- list(dual(3, 0), dual(4, 1))
  results <- Map("+", xs, ys)
  expect_equal(value(results[[1]]), 4)
  expect_equal(deriv(results[[1]]), 1)
  expect_equal(value(results[[2]]), 6)
  expect_equal(deriv(results[[2]]), 1)
})

test_that("vapply extracts values from list of duals", {
  duals <- list(dual(10, 1), dual(20, 2), dual(30, 3))
  vals <- vapply(duals, value, numeric(1))
  expect_equal(vals, c(10, 20, 30))
})

# -- Nested function composition -----------------------------------------------

test_that("nested user-defined functions compose correctly", {
  f <- function(x) x^2 + 1
  g <- function(x) sqrt(x)
  x <- dual_variable(3)
  # g(f(x)) = sqrt(x^2 + 1); d/dx = x / sqrt(x^2 + 1)
  r <- g(f(x))
  expect_equal(value(r), sqrt(10))
  expect_equal(deriv(r), 3 / sqrt(10), tolerance = 1e-12)
})

# -- while loop with dual accumulation ----------------------------------------

test_that("while loop accumulates dual correctly", {
  x <- dual_variable(2)
  accum <- dual(1, 0)
  i <- 0
  while (i < 3) {
    accum <- accum * x
    i <- i + 1
  }
  # accum = x^3 at x=2: value=8, deriv=3*x^2=12
  expect_equal(value(accum), 8)
  expect_equal(deriv(accum), 12)
})

# -- Dual in ifelse-like branches ---------------------------------------------

test_that("dual works with sequential if-else branching", {
  x <- dual(5, 1)
  result <- if (x > 3) {
    x^2
  } else if (x > 0) {
    x
  } else {
    -x
  }
  expect_equal(value(result), 25)
  expect_equal(deriv(result), 10)
})

# -- Mixed arithmetic chains --------------------------------------------------

test_that("long mixed arithmetic chain preserves derivatives", {
  x <- dual_variable(2)
  # ((x + 3) * 2 - 1) / x = (2x+5)/x = 2 + 5/x
  r <- ((x + 3) * 2 - 1) / x
  expect_equal(value(r), 9 / 2)
  # d/dx [2 + 5/x] = -5/x^2
  expect_equal(deriv(r), -5 / 4, tolerance = 1e-12)
})

# -- Using dual with seq-like index loop ---------------------------------------

test_that("lapply index loop with dual parameter", {
  mu <- dual_variable(0)
  data <- c(-1, 0, 1)
  # sum of squared residuals: sum((xi - mu)^2) at mu=0 = 1+0+1 = 2
  terms <- lapply(data, function(xi) (xi - mu)^2)
  ss <- Reduce("+", terms)
  expect_equal(value(ss), 2)
  # d/dmu sum((xi-mu)^2) = -2*sum(xi-mu) = -2*(-1+0+1) = 0
  expect_equal(deriv(ss), 0)
})

# -- dual_vector with sapply extraction ----------------------------------------

test_that("sapply on dual_vector elements", {
  dv <- dual_vector(dual(10, 1), dual(20, 2), dual(30, 3))
  vals <- sapply(seq_len(length(dv)), function(i) value(dv[i]))
  derivs <- sapply(seq_len(length(dv)), function(i) deriv(dv[i]))
  expect_equal(vals, c(10, 20, 30))
  expect_equal(derivs, c(1, 2, 3))
})

# -- Dual with do.call --------------------------------------------------------

test_that("do.call sum with dual args", {
  args <- list(dual(1, 1), dual(2, 0), dual(3, 0))
  r <- do.call(sum, args)
  expect_true(is_dual(r))
  expect_equal(value(r), 6)
  expect_equal(deriv(r), 1)
})

# -- Recursive function with dual ---------------------------------------------

test_that("recursive function with dual (Horner's method)", {
  # Evaluate polynomial 1 + 2x + 3x^2 at x using Horner's: ((3)*x + 2)*x + 1
  horner <- function(coeffs, x) {
    result <- dual(0, 0)
    for (co in rev(coeffs)) {
      result <- result * x + co
    }
    result
  }
  x <- dual_variable(2)
  r <- horner(c(1, 2, 3), x)
  # value: 1 + 4 + 12 = 17
  expect_equal(value(r), 17)
  # derivative: 2 + 6x = 2 + 12 = 14
  expect_equal(deriv(r), 14)
})
