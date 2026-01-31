test_that("dual constructor and accessors work", {
  d <- dual(3, 5)
  expect_equal(value(d), 3)
  expect_equal(deriv(d), 5)
})

test_that("dual_variable sets deriv=1", {
  d <- dual_variable(7)
  expect_equal(value(d), 7)
  expect_equal(deriv(d), 1)
})

test_that("dual_constant sets deriv=0", {
  d <- dual_constant(7)
  expect_equal(value(d), 7)
  expect_equal(deriv(d), 0)
})

test_that("is_dual predicate works", {
  expect_true(is_dual(dual(1, 0)))
  expect_false(is_dual(1))
  expect_false(is_dual("a"))
})

test_that("as.numeric extracts value", {
  expect_equal(as.numeric(dual(3.14, 1)), 3.14)
})

test_that("value/deriv on numeric pass through", {
  expect_equal(value(5), 5)
  expect_equal(deriv(5), 0)
})

# -- Addition ------------------------------------------------------------------

test_that("dual + dual", {
  a <- dual(3, 2)
  b <- dual(5, 7)
  r <- a + b
  expect_equal(value(r), 8)
  expect_equal(deriv(r), 9)
})

test_that("dual + numeric", {
  a <- dual(3, 2)
  r <- a + 10
  expect_equal(value(r), 13)
  expect_equal(deriv(r), 2)
})

test_that("numeric + dual", {
  a <- dual(3, 2)
  r <- 10 + a
  expect_equal(value(r), 13)
  expect_equal(deriv(r), 2)
})

# -- Subtraction ---------------------------------------------------------------

test_that("dual - dual", {
  a <- dual(10, 3)
  b <- dual(4, 1)
  r <- a - b
  expect_equal(value(r), 6)
  expect_equal(deriv(r), 2)
})

test_that("dual - numeric", {
  r <- dual(10, 3) - 4
  expect_equal(value(r), 6)
  expect_equal(deriv(r), 3)
})

test_that("numeric - dual", {
  r <- 10 - dual(4, 1)
  expect_equal(value(r), 6)
  expect_equal(deriv(r), -1)
})

# -- Multiplication ------------------------------------------------------------

test_that("dual * dual (product rule)", {
  # d/dx [x * x] at x=3 should be 2*3 = 6
  x <- dual(3, 1)
  r <- x * x
  expect_equal(value(r), 9)
  expect_equal(deriv(r), 6)
})

test_that("dual * numeric", {
  x <- dual(3, 1)
  r <- x * 5
  expect_equal(value(r), 15)
  expect_equal(deriv(r), 5)
})

test_that("numeric * dual", {
  x <- dual(3, 1)
  r <- 5 * x
  expect_equal(value(r), 15)
  expect_equal(deriv(r), 5)
})

test_that("0 * dual gives zero derivative", {
  x <- dual(3, 1)
  r <- 0 * x
  expect_equal(value(r), 0)
  expect_equal(deriv(r), 0)
})

# -- Division ------------------------------------------------------------------

test_that("dual / dual (quotient rule)", {
  # d/dx [x / (x+1)] = 1/(x+1)^2
  # At x=2: value = 2/3, deriv = 1/9
  x <- dual(2, 1)
  r <- x / (x + 1)
  expect_equal(value(r), 2/3)
  expect_equal(deriv(r), 1/9, tolerance = 1e-12)
})

test_that("dual / numeric", {
  x <- dual(6, 2)
  r <- x / 3
  expect_equal(value(r), 2)
  expect_equal(deriv(r), 2/3)
})

test_that("numeric / dual", {
  # d/dx [1/x] = -1/x^2
  # At x=2: value = 0.5, deriv = -0.25
  x <- dual(2, 1)
  r <- 1 / x
  expect_equal(value(r), 0.5)
  expect_equal(deriv(r), -0.25)
})

test_that("division by zero dual produces Inf", {
  r <- dual(1, 1) / dual(0, 1)
  expect_true(is.infinite(value(r)))
})

test_that("zero divided by zero dual produces NaN", {
  r <- dual(0, 1) / dual(0, 1)
  expect_true(is.nan(value(r)))
})

# -- Power ---------------------------------------------------------------------

test_that("dual ^ numeric (power rule)", {
  # d/dx [x^3] = 3x^2
  # At x=2: value = 8, deriv = 12
  x <- dual(2, 1)
  r <- x^3
  expect_equal(value(r), 8)
  expect_equal(deriv(r), 12)
})

test_that("numeric ^ dual (exponential rule)", {
  # d/dx [2^x] = 2^x * log(2)
  # At x=3: value = 8, deriv = 8*log(2)
  x <- dual(3, 1)
  r <- 2^x
  expect_equal(value(r), 8)
  expect_equal(deriv(r), 8 * log(2))
})

test_that("dual ^ dual (general power)", {
  # d/dx [x^x] = x^x * (log(x) + 1)
  # At x=2: value = 4, deriv = 4*(log(2)+1)
  x <- dual(2, 1)
  r <- x^x
  expect_equal(value(r), 4)
  expect_equal(deriv(r), 4 * (log(2) + 1))
})

# -- Unary operators -----------------------------------------------------------

test_that("unary minus", {
  x <- dual(3, 2)
  r <- -x
  expect_equal(value(r), -3)
  expect_equal(deriv(r), -2)
})

test_that("unary plus", {
  x <- dual(3, 2)
  r <- +x
  expect_equal(value(r), 3)
  expect_equal(deriv(r), 2)
})

# -- Comparison operators ------------------------------------------------------

test_that("comparison operators work on value", {
  a <- dual(3, 1)
  b <- dual(5, 2)
  expect_true(a < b)
  expect_true(a <= b)
  expect_false(a > b)
  expect_false(a >= b)
  expect_false(a == b)
  expect_true(a != b)
})

test_that("comparison with numeric", {
  x <- dual(3, 1)
  expect_true(x < 5)
  expect_true(x > 1)
  expect_true(x == 3)
  expect_true(2 < x)
})

# -- Chain rule verification ---------------------------------------------------

test_that("chain rule: d/dx [f(g(x))]", {
  # f(g(x)) = (2x+1)^3
  # f'(g(x)) * g'(x) = 3*(2x+1)^2 * 2
  # At x=1: value = 27, deriv = 3*9*2 = 54
  x <- dual(1, 1)
  r <- (2 * x + 1)^3
  expect_equal(value(r), 27)
  expect_equal(deriv(r), 54)
})

test_that("product and chain rules combined", {
  # d/dx [x^2 * (x+1)] = 2x*(x+1) + x^2 = 3x^2 + 2x
  # At x=3: value = 9*4 = 36, deriv = 27 + 6 = 33
  x <- dual(3, 1)
  r <- x^2 * (x + 1)
  expect_equal(value(r), 36)
  expect_equal(deriv(r), 33)
})

# -- dual_vector indexing ------------------------------------------------------

test_that("dual_vector supports [i] indexing", {
  dv <- dual_vector(dual(1, 0), dual(2, 1), dual(3, 0))
  expect_equal(length(dv), 3)
  d2 <- dv[2]
  expect_true(is_dual(d2))
  expect_equal(value(d2), 2)
  expect_equal(deriv(d2), 1)
})

test_that("dual_vector supports [[i]] indexing", {
  dv <- dual_vector(dual(10, 0), dual(20, 1))
  d <- dv[[1]]
  expect_true(is_dual(d))
  expect_equal(value(d), 10)
})
