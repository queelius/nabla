test_that("constant zero powers have zero derivatives at the origin", {
  expect_equal(D(function(x) x[1]^0, 0), 0)
  expect_equal(D(function(x) x[1]^dual_constant(0), 0), 0)
  expect_equal(deriv(0^dual_constant(2)), 0)
})

test_that("polynomial derivatives terminate at zero without NaN", {
  expect_equal(as.vector(D(function(x) x[1]^2, 0, order = 2)), 2)
  expect_equal(as.vector(D(function(x) x[1]^2, 0, order = 3)), 0)
  expect_equal(as.vector(D(function(x) x[1]^3, 0, order = 4)), 0)
  expect_equal(as.vector(D(function(x) x[1]^4, 0, order = 4)), 24)
})

test_that("zero first derivatives retain higher-order exponent curvature", {
  expect_equal(as.vector(D(function(x) 2^(x[1]^2), 0, order = 2)),
               2 * log(2))
  expect_equal(as.vector(D(function(x) (x[1] + 2)^(x[1]^2),
                          0, order = 2)), 2 * log(2))
})

test_that("scalar cumulative extrema preserve derivatives", {
  expect_equal(D(function(x) cummax(x[1]^2), 2), 4)
  expect_equal(D(function(x) cummin(x[1]^2), 2), 4)
  expect_equal(as.vector(D(function(x) cummax(x[1]^2), 2, order = 2)), 2)
  expect_equal(as.vector(D(function(x) cummin(x[1]^2), 2, order = 2)), 2)
  expect_error(cummax(dual(c(1, 2), c(1, 1))), "scalar")
  expect_error(cummin(dual(c(1, 2), c(1, 1))), "scalar")
  expect_error(cummax(dual(dual(c(1, 2), 0), 1)), "scalar")
})
