# Tests verifying all code examples from the three vignettes.
# Organized by vignette and section.

# ============================================================================
# Vignette 1: Introduction
# ============================================================================

test_that("basic dual creation and access", {
  x <- dual_variable(3)
  expect_equal(value(x), 3)
  expect_equal(deriv(x), 1)

  k <- dual_constant(5)
  expect_equal(value(k), 5)
  expect_equal(deriv(k), 0)

  y <- dual(2, 1)
  expect_equal(value(y), 2)
  expect_equal(deriv(y), 1)
})

test_that("arithmetic derivative rules", {
  x <- dual_variable(3)

  # Addition: d/dx(x + 2) = 1
  r <- x + 2
  expect_equal(value(r), 5)
  expect_equal(deriv(r), 1)

  # Subtraction: d/dx(5 - x) = -1
  r <- 5 - x
  expect_equal(value(r), 2)
  expect_equal(deriv(r), -1)

  # Multiplication: d/dx(x * 4) = 4
  r <- x * 4
  expect_equal(value(r), 12)
  expect_equal(deriv(r), 4)

  # Division: d/dx(1/x) = -1/x^2
  r <- 1 / x
  expect_equal(value(r), 1 / 3)
  expect_equal(deriv(r), -1 / 9)

  # Power: d/dx(x^3) = 3*x^2 = 27
  r <- x^3
  expect_equal(value(r), 27)
  expect_equal(deriv(r), 27)
})

test_that("math function derivatives", {
  # exp: d/dx exp(x) = exp(x)
  x <- dual_variable(1)
  r <- exp(x)
  expect_equal(value(r), exp(1))
  expect_equal(deriv(r), exp(1))

  # log: d/dx log(x) = 1/x
  r <- log(x)
  expect_equal(value(r), 0)
  expect_equal(deriv(r), 1)

  # sin: d/dx sin(x) = cos(x)
  x2 <- dual_variable(pi / 4)
  r <- sin(x2)
  expect_equal(value(r), sin(pi / 4))
  expect_equal(deriv(r), cos(pi / 4))

  # sqrt: d/dx sqrt(x) = 1/(2*sqrt(x))
  x3 <- dual_variable(4)
  r <- sqrt(x3)
  expect_equal(value(r), 2)
  expect_equal(deriv(r), 0.25)

  # lgamma: d/dx lgamma(x) = digamma(x)
  x4 <- dual_variable(3)
  r <- lgamma(x4)
  expect_equal(value(r), lgamma(3))
  expect_equal(deriv(r), digamma(3))
})

test_that("composition with user-defined functions", {
  f <- function(x) x^2 + sin(x)
  x <- dual_variable(pi / 4)
  result <- f(x)

  expected_val <- (pi / 4)^2 + sin(pi / 4)
  expected_deriv <- 2 * (pi / 4) + cos(pi / 4)

  expect_equal(value(result), expected_val)
  expect_equal(deriv(result), expected_deriv, tolerance = 1e-14)
})

test_that("standard normal PDF composition", {
  g <- function(x) exp(-x^2 / 2) / sqrt(2 * pi)
  x <- dual_variable(1)
  result <- g(x)

  expect_equal(value(result), dnorm(1), tolerance = 1e-14)
  expect_equal(deriv(result), -1 * dnorm(1), tolerance = 1e-14)
})

test_that("base R interop: sum, prod, c, is.numeric", {
  a <- dual_variable(2)
  b <- dual_constant(3)

  # sum
  total <- sum(a, b, dual_constant(1))
  expect_equal(value(total), 6)
  expect_equal(deriv(total), 1)

  # prod
  p <- prod(a, dual_constant(3))
  expect_equal(value(p), 6)
  expect_equal(deriv(p), 3)

  # c creates dual_vector
  v <- c(a, b)
  expect_s4_class(v, "dual_vector")
  expect_equal(length(v), 2)

  # is.numeric
  expect_true(is.numeric(dual_variable(1)))
})

test_that("if/else branching with duals", {
  safe_log <- function(x) {
    if (x > 0) log(x) else dual_constant(-Inf)
  }
  expect_equal(value(safe_log(dual_variable(2))), log(2))
  expect_equal(deriv(safe_log(dual_variable(2))), 0.5)
})

test_that("for loop accumulation", {
  x <- dual_variable(2)
  accum <- dual_constant(0)
  for (i in 1:5) {
    accum <- accum + x^i
  }
  # sum = 2 + 4 + 8 + 16 + 32 = 62
  expect_equal(value(accum), 62)
  # derivative of sum(x^i) = sum(i*x^(i-1)) = 1 + 4 + 12 + 32 + 80 = 129
  # at x=2: 1*2^0 + 2*2^1 + 3*2^2 + 4*2^3 + 5*2^4 = 1+4+12+32+80 = 129
  expect_equal(deriv(accum), 129)
})

test_that("Reduce with duals", {
  terms <- list(dual_variable(2), dual_constant(3), dual_constant(4))
  total <- Reduce("+", terms)
  expect_equal(value(total), 9)
  expect_equal(deriv(total), 1)
})

test_that("three-way comparison: x^3*sin(x)", {
  f <- function(x) x^3 * sin(x)
  x0 <- 2

  # Analytical
  analytical <- 3 * x0^2 * sin(x0) + x0^3 * cos(x0)

  # Finite differences
  h <- 1e-8
  finite_diff <- (f(x0 + h) - f(x0 - h)) / (2 * h)

  # AD
  ad_result <- deriv(f(dual_variable(x0)))

  # AD should match analytical to machine precision

  expect_equal(ad_result, analytical, tolerance = 1e-14)
  # Finite diff should be close but not exact
  expect_equal(finite_diff, analytical, tolerance = 1e-6)
})

# ============================================================================
# Vignette 2: MLE Workflow
# ============================================================================

# -- Normal(mu), known sigma --

test_that("Normal(mu) score equals analytical", {
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
  ad_s <- score(ll, mu0)
  analytical_s <- (sum_x - n * mu0) / sigma^2
  expect_equal(ad_s[1], analytical_s, tolerance = 1e-12)
})

test_that("Normal(mu) Hessian equals analytical", {
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
  ad_h <- hessian(ll, mu0)
  analytical_h <- -n / sigma^2
  expect_equal(ad_h[1, 1], analytical_h, tolerance = 1e-12)
})

test_that("Normal(mu) three-way agreement", {
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
  ad_s <- score(ll, mu0)
  ad_h <- hessian(ll, mu0)
  analytical_s <- (sum_x - n * mu0) / sigma^2
  analytical_h <- -n / sigma^2
  num_s <- numerical_gradient(ll, mu0)
  num_h <- numerical_hessian(ll, mu0)

  expect_equal(ad_s[1], analytical_s, tolerance = 1e-12)
  expect_equal(ad_h[1, 1], analytical_h, tolerance = 1e-12)
  expect_equal(num_s[1], analytical_s, tolerance = 1e-6)
  expect_equal(num_h[1, 1], analytical_h, tolerance = 1e-4)
})

# -- Normal(mu, sigma) --

test_that("Normal(mu, sigma) score vector matches analytical", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2
  xbar <- fix$xbar

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  theta0 <- c(4.5, 1.8)
  mu0 <- theta0[1]; sigma0 <- theta0[2]
  ss <- sum_x2 - 2 * mu0 * sum_x + n * mu0^2

  analytical_score <- c(
    n * (xbar - mu0) / sigma0^2,
    -n / sigma0 + ss / sigma0^3
  )

  ad_s <- score(ll, theta0)
  expect_equal(ad_s, analytical_score, tolerance = 1e-10)
})

test_that("Normal(mu, sigma) Hessian matrix matches analytical", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2
  xbar <- fix$xbar

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  theta0 <- c(4.5, 1.8)
  mu0 <- theta0[1]; sigma0 <- theta0[2]
  ss <- sum_x2 - 2 * mu0 * sum_x + n * mu0^2

  analytical_hess <- matrix(c(
    -n / sigma0^2,
    -2 * n * (xbar - mu0) / sigma0^3,
    -2 * n * (xbar - mu0) / sigma0^3,
    n / sigma0^2 - 3 * ss / sigma0^4
  ), nrow = 2, byrow = TRUE)

  ad_h <- hessian(ll, theta0)
  expect_equal(ad_h, analytical_hess, tolerance = 1e-10)
})

test_that("Normal(mu, sigma) three-way agreement", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  theta0 <- c(4.5, 1.8)

  ad_s <- score(ll, theta0)
  ad_h <- hessian(ll, theta0)
  num_s <- numerical_gradient(ll, theta0)
  num_h <- numerical_hessian(ll, theta0)

  expect_equal(ad_s, num_s, tolerance = 1e-6)
  expect_equal(ad_h, num_h, tolerance = 1e-4)
})

# -- Poisson(lambda) --

test_that("Poisson(lambda) score and Hessian match analytical", {
  set.seed(123)
  data_pois <- rpois(80, lambda = 3.5)
  n_pois <- length(data_pois)
  sum_x_pois <- sum(data_pois)
  sum_lfact <- sum(lfactorial(data_pois))

  ll <- function(theta) {
    lambda <- theta[1]
    sum_x_pois * log(lambda) - n_pois * lambda - sum_lfact
  }

  lam0 <- 3.0
  ad_s <- score(ll, lam0)
  ad_h <- hessian(ll, lam0)

  analytical_s <- sum_x_pois / lam0 - n_pois
  analytical_h <- -sum_x_pois / lam0^2

  expect_equal(ad_s[1], analytical_s, tolerance = 1e-12)
  expect_equal(ad_h[1, 1], analytical_h, tolerance = 1e-12)
})

test_that("Poisson(lambda) three-way agreement", {
  set.seed(123)
  data_pois <- rpois(80, lambda = 3.5)
  n_pois <- length(data_pois)
  sum_x_pois <- sum(data_pois)
  sum_lfact <- sum(lfactorial(data_pois))

  ll <- function(theta) {
    lambda <- theta[1]
    sum_x_pois * log(lambda) - n_pois * lambda - sum_lfact
  }

  lam0 <- 3.0
  ad_s <- score(ll, lam0)
  ad_h <- hessian(ll, lam0)
  num_s <- numerical_gradient(ll, lam0)
  num_h <- numerical_hessian(ll, lam0)
  analytical_s <- sum_x_pois / lam0 - n_pois
  analytical_h <- -sum_x_pois / lam0^2

  expect_equal(ad_s[1], analytical_s, tolerance = 1e-12)
  expect_equal(num_s[1], analytical_s, tolerance = 1e-6)
  expect_equal(ad_h[1, 1], analytical_h, tolerance = 1e-12)
  expect_equal(num_h[1, 1], analytical_h, tolerance = 1e-4)
})

# -- Gamma(shape), known rate --

test_that("Gamma(shape) score and Hessian match analytical", {
  set.seed(99)
  data_gamma <- rgamma(60, shape = 2.5, rate = 1)
  n_gam <- length(data_gamma)
  sum_log_x <- sum(log(data_gamma))
  sum_x_gam <- sum(data_gamma)
  beta_known <- 1

  ll <- function(theta) {
    alpha <- theta[1]
    (alpha - 1) * sum_log_x - n_gam * lgamma(alpha) +
      n_gam * alpha * log(beta_known) - beta_known * sum_x_gam
  }

  alpha0 <- 2.0
  ad_s <- score(ll, alpha0)
  ad_h <- hessian(ll, alpha0)

  analytical_s <- sum_log_x - n_gam * digamma(alpha0) + n_gam * log(beta_known)
  analytical_h <- -n_gam * trigamma(alpha0)

  expect_equal(ad_s[1], analytical_s, tolerance = 1e-10)
  expect_equal(ad_h[1, 1], analytical_h, tolerance = 1e-10)
})

test_that("Gamma(shape) three-way agreement", {
  set.seed(99)
  data_gamma <- rgamma(60, shape = 2.5, rate = 1)
  n_gam <- length(data_gamma)
  sum_log_x <- sum(log(data_gamma))
  sum_x_gam <- sum(data_gamma)
  beta_known <- 1

  ll <- function(theta) {
    alpha <- theta[1]
    (alpha - 1) * sum_log_x - n_gam * lgamma(alpha) +
      n_gam * alpha * log(beta_known) - beta_known * sum_x_gam
  }

  alpha0 <- 2.0
  ad_s <- score(ll, alpha0)
  ad_h <- hessian(ll, alpha0)
  num_s <- numerical_gradient(ll, alpha0)
  num_h <- numerical_hessian(ll, alpha0)

  expect_equal(ad_s[1], num_s[1], tolerance = 1e-6)
  expect_equal(ad_h[1, 1], num_h[1, 1], tolerance = 1e-4)
})

# -- Logistic regression --

test_that("Logistic regression score and Hessian match numerical", {
  set.seed(7)
  n_lr <- 50
  X <- cbind(1, rnorm(n_lr), rnorm(n_lr))
  beta_true <- c(-0.5, 1.2, -0.8)
  eta_true <- X %*% beta_true
  prob_true <- 1 / (1 + exp(-eta_true))
  y <- rbinom(n_lr, 1, prob_true)

  ll_dual <- function(theta) {
    result <- dual_constant(0)
    for (i in seq_len(n_lr)) {
      eta_i <- theta[1] * X[i, 1] + theta[2] * X[i, 2] + theta[3] * X[i, 3]
      result <- result + y[i] * eta_i - log(1 + exp(eta_i))
    }
    result
  }

  ll_num <- function(beta) {
    eta <- X %*% beta
    sum(y * eta - log(1 + exp(eta)))
  }

  beta0 <- c(0, 0, 0)
  ad_s <- score(ll_dual, beta0)
  ad_h <- hessian(ll_dual, beta0)
  num_s <- numerical_gradient(ll_num, beta0)
  num_h <- numerical_hessian(ll_num, beta0)

  expect_equal(ad_s, num_s, tolerance = 1e-6)
  expect_equal(ad_h, num_h, tolerance = 1e-4)
})

test_that("Logistic regression three-way agreement", {
  set.seed(7)
  n_lr <- 50
  X <- cbind(1, rnorm(n_lr), rnorm(n_lr))
  beta_true <- c(-0.5, 1.2, -0.8)
  eta_true <- X %*% beta_true
  prob_true <- 1 / (1 + exp(-eta_true))
  y <- rbinom(n_lr, 1, prob_true)

  ll_dual <- function(theta) {
    result <- dual_constant(0)
    for (i in seq_len(n_lr)) {
      eta_i <- theta[1] * X[i, 1] + theta[2] * X[i, 2] + theta[3] * X[i, 3]
      result <- result + y[i] * eta_i - log(1 + exp(eta_i))
    }
    result
  }

  ll_num <- function(beta) {
    eta <- X %*% beta
    sum(y * eta - log(1 + exp(eta)))
  }

  beta0 <- c(0, 0, 0)
  ad_s <- score(ll_dual, beta0)
  ad_h <- hessian(ll_dual, beta0)
  num_s <- numerical_gradient(ll_num, beta0)
  num_h <- numerical_hessian(ll_num, beta0)

  # Analytical for logistic at beta=0: score = X'(y - 0.5), Hessian = -0.25 * X'X
  p_hat <- rep(0.5, n_lr)
  analytical_s <- as.numeric(t(X) %*% (y - p_hat))
  analytical_h <- -t(X) %*% diag(p_hat * (1 - p_hat)) %*% X

  expect_equal(ad_s, analytical_s, tolerance = 1e-10)
  expect_equal(ad_h, unname(analytical_h), tolerance = 1e-10)
  expect_equal(num_s, analytical_s, tolerance = 1e-5)
  expect_equal(num_h, unname(analytical_h), tolerance = 1e-3)
})

# -- Newton-Raphson --

test_that("Newton-Raphson converges to MLE for Normal model", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  result <- newton_raphson(ll, c(3, 1))

  expect_equal(result$estimate[1], fix$mle_mu, tolerance = 1e-6)
  expect_equal(result$estimate[2], fix$mle_sigma, tolerance = 1e-6)
})

# -- score_and_hessian() --

test_that("score_and_hessian() matches hessian() from log-likelihood", {
  fix <- make_normal_fixture()
  n <- fix$n
  sum_x <- fix$sum_x
  sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  score_fn <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    s_mu <- (sum_x - n * mu) / sigma^2
    s_sigma <- -n / sigma + (sum_x2 - 2 * mu * sum_x + n * mu^2) / sigma^3
    list(s_mu, s_sigma)
  }

  theta0 <- c(4.5, 1.8)

  result_sh <- score_and_hessian(score_fn, theta0)
  hess_from_ll <- hessian(ll, theta0)

  expect_equal(result_sh$hessian, hess_from_ll, tolerance = 1e-10)
  expect_equal(result_sh$score, score(ll, theta0), tolerance = 1e-10)
})

# ============================================================================
# Vignette 3: Higher-Order Derivatives
# ============================================================================

test_that("differentiate2() returns correct value, first, second derivatives", {
  # sin(x) at pi/4
  result <- differentiate2(sin, pi / 4)
  expect_equal(result$value, sin(pi / 4), tolerance = 1e-14)
  expect_equal(result$first, cos(pi / 4), tolerance = 1e-14)
  expect_equal(result$second, -sin(pi / 4), tolerance = 1e-14)
})

test_that("differentiate2() for x*exp(-x^2)", {
  f <- function(x) x * exp(-x^2)
  d2 <- differentiate2(f, 1)

  # f'(x) = exp(-x^2)(1 - 2x^2)
  # f''(x) = exp(-x^2)(-6x + 4x^3)
  expect_equal(d2$value, exp(-1), tolerance = 1e-14)
  expect_equal(d2$first, exp(-1) * (1 - 2), tolerance = 1e-14)
  expect_equal(d2$second, exp(-1) * (-6 + 4), tolerance = 1e-14)
})

test_that("curvature formula gives correct result", {
  curvature <- function(f, x) {
    d2 <- differentiate2(f, x)
    abs(d2$second) / (1 + d2$first^2)^(3 / 2)
  }

  # sin(x) at pi/2: f'=0, f''=-1, kappa = 1/(1+0)^(3/2) = 1
  expect_equal(curvature(sin, pi / 2), 1.0, tolerance = 1e-14)

  # sin(x) at 0: f'=1, f''=0, kappa = 0
  expect_equal(curvature(sin, 0), 0, tolerance = 1e-14)
})

test_that("Taylor approximation is accurate near expansion point", {
  taylor2 <- function(f, x0, x) {
    d2 <- differentiate2(f, x0)
    d2$value + d2$first * (x - x0) + 0.5 * d2$second * (x - x0)^2
  }

  # exp(x) around x=0
  expect_equal(taylor2(exp, 0, 0), 1)
  expect_equal(taylor2(exp, 0, 0.01), exp(0.01), tolerance = 1e-6)
  expect_equal(taylor2(exp, 0, -0.01), exp(-0.01), tolerance = 1e-6)

  # Larger displacement: should still be close
  expect_equal(taylor2(exp, 0, 0.1), exp(0.1), tolerance = 1e-3)
})

test_that("manual nested dual approach matches hessian() helper", {
  set.seed(123)
  data_pois <- rpois(50, lambda = 3)
  n <- length(data_pois)
  sum_x <- sum(data_pois)
  sum_lfact <- sum(lfactorial(data_pois))

  ll <- function(theta) {
    lambda <- theta[1]
    sum_x * log(lambda) - n * lambda - sum_lfact
  }

  lambda0 <- 2.5

  # hessian() helper
  hess_helper <- hessian(ll, lambda0)

  # manual nested dual
  manual_theta <- dual_vector(list(dual2_variable(lambda0)))
  result_manual <- ll(manual_theta)
  manual_hess <- deriv(deriv(result_manual))

  expect_equal(hess_helper[1, 1], manual_hess, tolerance = 1e-14)
})

test_that("dual2_variable basic extraction", {
  x <- dual2_variable(2)
  result <- x^3

  expect_equal(value2(result), 8)
  expect_equal(first_deriv(result), 12)
  expect_equal(second_deriv(result), 12)
})

test_that("dual2_constant has zero derivatives", {
  k <- dual2_constant(5)
  expect_equal(value2(k), 5)
  expect_equal(first_deriv(k), 0)
  expect_equal(second_deriv(k), 0)
})

# ============================================================================
# Plot data generation tests
# ============================================================================

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

test_that("mle plot 1: Poisson LL + score grid produces finite values", {
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
  sc_vals <- sapply(lam_grid, function(l) score(ll, l))

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
    d2 <- differentiate2(f, x0)
    d2$value + d2$first * (x - x0) + 0.5 * d2$second * (x - x0)^2
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
    d2 <- differentiate2(f, x)
    abs(d2$second) / (1 + d2$first^2)^(3 / 2)
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
  ngr <- function(theta) -score(ll, theta)

  # BFGS path collection
  trace <- list()
  fn_trace <- function(theta) {
    trace[[length(trace) + 1L]] <<- theta
    nll(theta)
  }
  trace <- list()
  optim(c(0, 1), fn = fn_trace, gr = ngr, method = "BFGS")
  path <- do.call(rbind, trace)

  expect_true(all(is.finite(path)))
  expect_true(nrow(path) >= 2)
})

test_that("observed_information returns negative Hessian", {
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
  obs_info <- observed_information(ll, mu0)
  hess <- hessian(ll, mu0)

  expect_equal(obs_info, -hess)
  expect_equal(obs_info[1, 1], n / sigma^2)
})
