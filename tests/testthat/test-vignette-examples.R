# Tests for vignette-specific workflows not covered by unit tests.
# Plot data generation, optimizer paths, contour surfaces, and
# observed information — these validate the vignette plotting code.

test_that("introduction plot 1: function + derivative grid produces finite values", {
  f <- function(x) x^2 + sin(x)
  xs <- seq(0, 2 * pi, length.out = 50)

  vals <- sapply(xs, function(xi) {
    r <- f(dual_variable(xi))
    c(value(r), deriv(r))
  })

  expect_true(all(is.finite(vals[1, ])))
  expect_true(all(is.finite(vals[2, ])))

  # Check specific point
  x_mark <- pi / 4
  r_mark <- f(dual_variable(x_mark))
  expect_equal(value(r_mark), (pi / 4)^2 + sin(pi / 4), tolerance = 1e-14)
  expect_equal(deriv(r_mark), 2 * (pi / 4) + cos(pi / 4), tolerance = 1e-14)
})

test_that("introduction plot 2: error comparison values are finite and ordered", {
  f <- function(x) x^3 * sin(x)
  x0 <- 2
  analytical <- 3 * x0^2 * sin(x0) + x0^3 * cos(x0)
  h <- 1e-8
  finite_diff <- (f(x0 + h) - f(x0 - h)) / (2 * h)
  ad_result <- deriv(f(dual_variable(x0)))

  err_fd <- abs(finite_diff - analytical)
  err_ad <- abs(ad_result - analytical)

  expect_true(is.finite(err_fd))
  expect_true(is.finite(err_ad))
  # AD error should be smaller than finite diff error
  expect_true(err_ad < err_fd)
})

test_that("mle plot 1: Poisson LL + gradient grid produces finite values", {
  set.seed(123)
  data_pois <- rpois(80, lambda = 3.5)
  n_pois <- length(data_pois)
  sum_x_pois <- sum(data_pois)
  sum_lfact <- sum(lfactorial(data_pois))

  ll <- function(theta) {
    lambda <- theta[1]
    sum_x_pois * log(lambda) - n_pois * lambda - sum_lfact
  }

  lam_grid <- seq(2.0, 5.5, length.out = 50)
  ll_vals <- sapply(lam_grid, function(l) ll(l))
  sc_vals <- sapply(lam_grid, function(l) gradient(ll, l))

  expect_true(all(is.finite(ll_vals)))
  expect_true(all(is.finite(sc_vals)))
})

test_that("mle plot 2: Normal contour surface produces finite values", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  mu_grid <- seq(4.0, 6.0, length.out = 10)
  sigma_grid <- seq(1.2, 2.8, length.out = 10)
  ll_surface <- outer(mu_grid, sigma_grid, Vectorize(function(m, s) {
    ll(c(m, s))
  }))

  expect_true(all(is.finite(ll_surface)))
})

test_that("mle plot 3: Newton-Raphson trace records valid iterates", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  result <- newton_raphson(ll, c(3, 1), trace = TRUE)

  expect_true(all(is.finite(result$trace)))
  expect_true(nrow(result$trace) >= 2)  # at least start + 1 step
  expect_true(result$iterations <= 50)

  # Final iterate should be near MLE
  expect_equal(result$estimate[1], fix$mle_mu, tolerance = 1e-6)
  expect_equal(result$estimate[2], fix$mle_sigma, tolerance = 1e-6)
})

test_that("higher-order plot 1: Taylor vs exact grid is finite", {
  taylor2 <- function(f, x0, x) {
    d2 <- differentiate_n(f, x0, order = 2)
    d2$value + d2$d1 * (x - x0) + 0.5 * d2$d2 * (x - x0)^2
  }

  xs <- seq(-2, 3, length.out = 50)
  exact_vals <- exp(xs)
  taylor_vals <- sapply(xs, function(x) taylor2(exp, 0, x))

  expect_true(all(is.finite(exact_vals)))
  expect_true(all(is.finite(taylor_vals)))

  # Near x0=0, Taylor should be very close
  expect_equal(taylor2(exp, 0, 0.01), exp(0.01), tolerance = 1e-6)
})

test_that("higher-order plot 2: sin curvature grid is finite", {
  curvature <- function(f, x) {
    d2 <- differentiate_n(f, x, order = 2)
    abs(d2$d2) / (1 + d2$d1^2)^(3 / 2)
  }

  xs <- seq(0, 2 * pi, length.out = 50)
  sin_vals <- sin(xs)
  kappa_vals <- sapply(xs, function(x) curvature(sin, x))

  expect_true(all(is.finite(sin_vals)))
  expect_true(all(is.finite(kappa_vals)))
  expect_true(all(kappa_vals >= 0))

  # Max curvature at pi/2
  expect_equal(curvature(sin, pi / 2), 1.0, tolerance = 1e-14)
})

test_that("optimizer plot: contour + path data is finite", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  nll <- function(theta) -ll(theta)
  ngr <- function(theta) -gradient(ll, theta)

  # BFGS path collection (use environment to avoid <<-)
  env <- new.env(parent = emptyenv())
  env$trace <- list()
  fn_trace <- function(theta) {
    env$trace[[length(env$trace) + 1L]] <- theta
    nll(theta)
  }
  optim(c(0, 1), fn = fn_trace, gr = ngr, method = "BFGS")
  path <- do.call(rbind, env$trace)

  expect_true(all(is.finite(path)))
  expect_true(nrow(path) >= 2)
})

test_that("-hessian() gives observed information", {
  fix <- make_normal_fixture()
  sigma <- 2
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]
    -1 / (2 * sigma^2) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  mu0 <- 4.5
  obs_info <- -hessian(ll, mu0)
  hess <- hessian(ll, mu0)

  expect_equal(obs_info, -hess)
  expect_equal(obs_info[1, 1], n / sigma^2)
})
