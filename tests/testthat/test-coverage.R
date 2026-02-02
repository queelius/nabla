# Tests for uncovered code paths across nabla
# Targets: modulo, integer division, logical ops, summary group,
# math functions, class display, special functions, standalone math.

# ===========================================================================
# 1. Arithmetic: modulo, integer division, logical operators
# ===========================================================================

# -- Modulo (%%) ---------------------------------------------------------------

test_that("%% (modulo) on dual,dual", {
  r <- dual(7, 1) %% dual(3, 0)
  expect_true(is_dual(r))
  expect_equal(value(r), 1)
  # derivative follows left operand's derivative

  expect_equal(deriv(r), 1)
})

test_that("%% (modulo) on dual,numeric", {
  r <- dual(7, 1) %% 3
  expect_true(is_dual(r))
  expect_equal(value(r), 1)
  expect_equal(deriv(r), 1)
})

test_that("%% (modulo) on numeric,dual", {
  r <- 7 %% dual(3, 1)
  expect_true(is_dual(r))
  expect_equal(value(r), 1)
})

# -- Integer division (%/%) ---------------------------------------------------

test_that("%/% (integer division) on dual,dual", {
  r <- dual(7, 1) %/% dual(3, 0)
  expect_true(is_dual(r))
  expect_equal(value(r), 2)
  # integer division is piecewise constant, derivative is 0
  expect_equal(deriv(r), 0)
})

test_that("%/% (integer division) on dual,numeric", {
  r <- dual(7, 1) %/% 3
  expect_true(is_dual(r))
  expect_equal(value(r), 2)
  expect_equal(deriv(r), 0)
})

test_that("%/% (integer division) on numeric,dual", {
  r <- 7 %/% dual(3, 1)
  expect_true(is_dual(r))
  expect_equal(value(r), 2)
  expect_equal(deriv(r), 0)
})

# -- Logical AND (&) ----------------------------------------------------------

test_that("& (logical AND) on dual,dual", {
  expect_true(dual(1, 1) & dual(1, 0))
  expect_false(dual(0, 1) & dual(1, 0))
  expect_false(dual(1, 1) & dual(0, 0))
  expect_false(dual(0, 0) & dual(0, 1))
})

test_that("& (logical AND) on dual,numeric", {
  expect_true(dual(1, 1) & 1)
  expect_false(dual(0, 1) & 1)
})

test_that("& (logical AND) on numeric,dual", {
  expect_true(1 & dual(1, 0))
  expect_false(0 & dual(1, 0))
})

# -- Logical OR (|) -----------------------------------------------------------

test_that("| (logical OR) on dual,dual", {
  expect_true(dual(0, 1) | dual(1, 0))
  expect_true(dual(1, 1) | dual(0, 0))
  expect_true(dual(1, 1) | dual(1, 0))
  expect_false(dual(0, 1) | dual(0, 0))
})

test_that("| (logical OR) on dual,numeric", {
  expect_true(dual(0, 1) | 1)
  expect_false(dual(0, 1) | 0)
})

test_that("| (logical OR) on numeric,dual", {
  expect_true(0 | dual(1, 0))
  expect_false(0 | dual(0, 0))
})

# -- Logical NOT (!) ----------------------------------------------------------

test_that("! (logical NOT) on dual", {
  expect_true(!dual(0, 1))
  expect_false(!dual(1, 0))
  expect_false(!dual(5, 3))
  expect_true(!dual(0, 0))
})

# ===========================================================================
# 2. Summary group: range, any, all
# ===========================================================================

# -- range() -------------------------------------------------------------------

test_that("range() on duals returns dual_vector with min and max", {
  r <- range(dual(3, 1), dual(1, 2))
  expect_true(is(r, "dual_vector"))
  expect_equal(length(r), 2)
  # First element is the min
  expect_equal(value(r[1]), 1)
  expect_equal(deriv(r[1]), 2)
  # Second element is the max
  expect_equal(value(r[2]), 3)
  expect_equal(deriv(r[2]), 1)
})

test_that("range() on three duals", {
  r <- range(dual(5, 0), dual(1, 1), dual(3, 2))
  expect_equal(value(r[1]), 1)
  expect_equal(deriv(r[1]), 1)
  expect_equal(value(r[2]), 5)
  expect_equal(deriv(r[2]), 0)
})

test_that("range() on dual and numeric", {
  r <- range(dual(3, 1), 1)
  expect_true(is(r, "dual_vector"))
  expect_equal(value(r[1]), 1)
  expect_equal(deriv(r[1]), 0)
  expect_equal(value(r[2]), 3)
  expect_equal(deriv(r[2]), 1)
})

# -- any() ---------------------------------------------------------------------

test_that("any() on duals returns logical", {
  expect_true(any(dual(0, 1), dual(1, 0)))
  expect_true(any(dual(1, 1), dual(1, 0)))
  expect_false(any(dual(0, 1), dual(0, 0)))
})

test_that("any() single TRUE dual", {
  expect_true(any(dual(1, 0)))
})

test_that("any() single FALSE dual", {
  expect_false(any(dual(0, 1)))
})

# -- all() ---------------------------------------------------------------------

test_that("all() on duals returns logical", {
  expect_true(all(dual(1, 1), dual(1, 0)))
  expect_false(all(dual(0, 1), dual(1, 0)))
  expect_false(all(dual(1, 1), dual(0, 0)))
})

test_that("all() single TRUE dual", {
  expect_true(all(dual(1, 0)))
})

test_that("all() single FALSE dual", {
  expect_false(all(dual(0, 1)))
})

# ===========================================================================
# 3. Math functions: trunc, round, cospi, sinpi, tanpi, cumsum,
#    cummax, cummin, factorial, lfactorial
# ===========================================================================

# -- trunc() -------------------------------------------------------------------

test_that("trunc() on dual", {
  r <- trunc(dual(3.7, 1))
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 0)
})

test_that("trunc() on negative dual", {
  r <- trunc(dual(-3.7, 1))
  expect_equal(value(r), -3)
  expect_equal(deriv(r), 0)
})

# -- round() -------------------------------------------------------------------

test_that("round() on dual", {
  r <- round(dual(3.7, 1))
  expect_true(is_dual(r))
  expect_equal(value(r), 4)
  expect_equal(deriv(r), 0)
})

test_that("round() on dual with value below .5", {
  r <- round(dual(3.2, 1))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 0)
})

# -- cospi() -------------------------------------------------------------------

test_that("cospi() on dual", {
  r <- cospi(dual(0.5, 1))
  expect_true(is_dual(r))
  # cos(pi/2) is approximately 0
  expect_equal(value(r), cospi(0.5), tolerance = 1e-13)
  # derivative: -pi * sin(pi/2) * 1 = -pi
  expect_equal(deriv(r), -pi * sinpi(0.5), tolerance = 1e-12)
})

test_that("cospi() at zero", {
  r <- cospi(dual(0, 1))
  # cos(0) = 1
  expect_equal(value(r), 1)
  # -pi * sin(0) = 0
  expect_equal(deriv(r), 0, tolerance = 1e-13)
})

# -- sinpi() -------------------------------------------------------------------

test_that("sinpi() on dual", {
  r <- sinpi(dual(0.5, 1))
  expect_true(is_dual(r))
  # sin(pi/2) = 1
  expect_equal(value(r), sinpi(0.5), tolerance = 1e-13)
  # derivative: pi * cos(pi/2) * 1 approximately 0
  expect_equal(deriv(r), pi * cospi(0.5), tolerance = 1e-12)
})

test_that("sinpi() at zero", {
  r <- sinpi(dual(0, 1))
  # sin(0) = 0
  expect_equal(value(r), 0, tolerance = 1e-13)
  # pi * cos(0) = pi
  expect_equal(deriv(r), pi, tolerance = 1e-12)
})

# -- tanpi() -------------------------------------------------------------------

test_that("tanpi() on dual", {
  r <- tanpi(dual(0.25, 1))
  expect_true(is_dual(r))
  # tan(pi/4) = 1
  expect_equal(value(r), tanpi(0.25), tolerance = 1e-12)
  # derivative: pi / cos(pi/4)^2 = pi / 0.5 = 2*pi
  cv <- cospi(0.25)
  expect_equal(deriv(r), pi / (cv * cv), tolerance = 1e-12)
})

test_that("tanpi() at zero", {
  r <- tanpi(dual(0, 1))
  # tan(0) = 0
  expect_equal(value(r), 0, tolerance = 1e-13)
  # pi / cos(0)^2 = pi
  expect_equal(deriv(r), pi, tolerance = 1e-12)
})

# -- cumsum() ------------------------------------------------------------------

test_that("cumsum() on scalar dual", {
  r <- cumsum(dual(3, 1))
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 1)
})

# -- cummax() ------------------------------------------------------------------

test_that("cummax() on scalar dual", {
  r <- cummax(dual(3, 1))
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 0)
})

# -- cummin() ------------------------------------------------------------------

test_that("cummin() on scalar dual", {
  r <- cummin(dual(3, 1))
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 0)
})

# -- factorial() ---------------------------------------------------------------

test_that("factorial() on dual", {
  r <- factorial(dual(4, 1))
  expect_true(is_dual(r))
  # factorial(4) = gamma(5) = 24
  expect_equal(value(r), gamma(5))
  # derivative: gamma(5) * digamma(5) * 1 = 24 * digamma(5)
  expect_equal(deriv(r), 24 * digamma(5), tolerance = 1e-10)
})

test_that("factorial() on dual with non-integer value", {
  r <- factorial(dual(2.5, 1))
  expect_equal(value(r), gamma(3.5), tolerance = 1e-10)
  expect_equal(deriv(r), gamma(3.5) * digamma(3.5), tolerance = 1e-10)
})

# -- lfactorial() --------------------------------------------------------------

test_that("lfactorial() on dual", {
  r <- lfactorial(dual(4, 1))
  expect_true(is_dual(r))
  # lfactorial(4) = lgamma(5) = log(24)
  expect_equal(value(r), lfactorial(4))
  # derivative: digamma(5)
  expect_equal(deriv(r), digamma(5), tolerance = 1e-10)
})

test_that("lfactorial() on dual with small value", {
  r <- lfactorial(dual(1, 1))
  expect_equal(value(r), lfactorial(1))
  expect_equal(deriv(r), digamma(2), tolerance = 1e-10)
})

# ===========================================================================
# 4. Class: show() display and dual_vector multi-index
# ===========================================================================

# -- show() for dual -----------------------------------------------------------

test_that("show() for dual produces expected output", {
  x <- dual(3, 1)
  out <- capture.output(show(x))
  expect_length(out, 1)
  expect_match(out, "dual")
  expect_match(out, "3")
  expect_match(out, "1")
})

# -- show() for dual_vector ----------------------------------------------------

test_that("show() for dual_vector produces expected output", {
  dv <- dual_vector(dual(1, 0), dual(2, 1))
  out <- capture.output(show(dv))
  # First line: header with element count
  expect_match(out[1], "dual_vector")
  expect_match(out[1], "2")
  # Subsequent lines: individual elements
  expect_length(out, 3)  # 1 header + 2 elements
  expect_match(out[2], "\\[1\\]")
  expect_match(out[3], "\\[2\\]")
})

test_that("show() for empty dual_vector", {
  dv <- dual_vector(list())
  out <- capture.output(show(dv))
  expect_match(out[1], "0")
})

# -- dual_vector multi-index ---------------------------------------------------

test_that("dual_vector multi-index [1:2] returns dual_vector", {
  dv <- dual_vector(dual(10, 1), dual(20, 2), dual(30, 3))
  sub <- dv[1:2]
  expect_true(is(sub, "dual_vector"))
  expect_equal(length(sub), 2)
  expect_equal(value(sub[1]), 10)
  expect_equal(deriv(sub[1]), 1)
  expect_equal(value(sub[2]), 20)
  expect_equal(deriv(sub[2]), 2)
})

test_that("dual_vector multi-index [c(1,3)]", {
  dv <- dual_vector(dual(10, 1), dual(20, 2), dual(30, 3))
  sub <- dv[c(1, 3)]
  expect_true(is(sub, "dual_vector"))
  expect_equal(length(sub), 2)
  expect_equal(value(sub[1]), 10)
  expect_equal(value(sub[2]), 30)
})

test_that("dual_vector single-index [2] returns dual (not dual_vector)", {
  dv <- dual_vector(dual(10, 1), dual(20, 2), dual(30, 3))
  elem <- dv[2]
  expect_true(is_dual(elem))
  expect_false(is(elem, "dual_vector"))
  expect_equal(value(elem), 20)
  expect_equal(deriv(elem), 2)
})

# ===========================================================================
# 5. Special: beta(numeric, dual)
# ===========================================================================

test_that("beta(numeric, dual) returns dual with correct value", {
  r <- beta(2, dual(3, 1))
  expect_true(is_dual(r))
  expect_equal(value(r), base::beta(2, 3), tolerance = 1e-12)
})

test_that("beta(numeric, dual) derivative matches numerical", {
  f <- function(b) beta(2, b)
  b_val <- 3
  r <- beta(2, dual(b_val, 1))
  num <- central_difference(f, b_val)
  expect_equal(deriv(r), num, tolerance = 1e-6)
})

test_that("beta(numeric, dual) at non-integer values", {
  r <- beta(0.5, dual(0.5, 1))
  expect_equal(value(r), base::beta(0.5, 0.5), tolerance = 1e-12)
  f <- function(b) beta(0.5, b)
  num <- central_difference(f, 0.5)
  expect_equal(deriv(r), num, tolerance = 1e-6)
})

# ===========================================================================
# 6. Standalone math: atan2(numeric,dual), min/max single-arg and mixed
# ===========================================================================

# -- atan2(numeric, dual) ------------------------------------------------------

test_that("atan2(numeric, dual) returns correct value and derivative", {
  r <- atan2(1, dual(2, 1))
  expect_true(is_dual(r))
  expect_equal(value(r), atan2(1, 2))
  # d/dx atan2(1, x) = -1/(1 + x^2) ... actually:
  # atan2(y,x) derivative wrt x: -y/(x^2 + y^2)
  # At y=1, x=2: -1/(4+1) = -1/5
  expect_equal(deriv(r), -1 / 5, tolerance = 1e-10)
})

test_that("atan2(numeric, dual) numerical check", {
  f <- function(x) atan2(1, x)
  x_val <- 2
  r <- atan2(1, dual(x_val, 1))
  num <- central_difference(f, x_val)
  expect_equal(deriv(r), num, tolerance = 1e-6)
})

# -- min single-arg ------------------------------------------------------------

test_that("min(dual) single argument returns the dual itself", {
  x <- dual(3, 1)
  r <- min(x)
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 1)
})

# -- max single-arg ------------------------------------------------------------

test_that("max(dual) single argument returns the dual itself", {
  x <- dual(3, 1)
  r <- max(x)
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 1)
})

# -- min(dual, numeric) --------------------------------------------------------

test_that("min(dual, numeric) where dual is smaller", {
  r <- min(dual(3, 1), 5)
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 1)
})

test_that("min(dual, numeric) where numeric is smaller", {
  r <- min(dual(3, 1), 1)
  expect_true(is_dual(r))
  expect_equal(value(r), 1)
  expect_equal(deriv(r), 0)
})

# -- max(dual, numeric) --------------------------------------------------------

test_that("max(dual, numeric) where dual is larger", {
  r <- max(dual(3, 1), 1)
  expect_true(is_dual(r))
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 1)
})

test_that("max(dual, numeric) where numeric is larger", {
  r <- max(dual(3, 1), 5)
  expect_true(is_dual(r))
  expect_equal(value(r), 5)
  expect_equal(deriv(r), 0)
})

# -- min/max with two duals (Summary group path) -------------------------------

test_that("min via Summary group selects correct dual", {
  a <- dual(7, 1)
  b <- dual(2, 3)
  # Use sum() style call to go through Summary group
  r <- min(a, b)
  expect_equal(value(r), 2)
  expect_equal(deriv(r), 3)
})

test_that("max via Summary group selects correct dual", {
  a <- dual(7, 1)
  b <- dual(2, 3)
  r <- max(a, b)
  expect_equal(value(r), 7)
  expect_equal(deriv(r), 1)
})

