# Tests for MLE helper functions: score, hessian, observed_information,
# score_and_hessian.
#
# Each test defines a log-likelihood with known analytical derivatives and
# verifies the AD results.
#
# IMPORTANT: Log-likelihood functions must keep dual parameters as duals
# throughout. Since base R's sum() doesn't dispatch on dual objects,
# compute sums algebraically (e.g., n*mu^2 - 2*mu*sum(data)) rather
# than element-wise (sum((data - mu)^2)).

tol_mle <- 1e-6

# -- Normal distribution (mu only) --------------------------------------------

test_that("Normal (mu only): score and hessian", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  # Log-likelihood keeping mu as dual:
  # -n/2*log(2pi) - 0.5 * (sum(x^2) - 2*mu*sum(x) + n*mu^2)
  loglik <- function(theta) {
    mu <- theta[1]
    -n/2 * log(2 * pi) - 0.5 * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  theta <- c(2)

  s <- score(loglik, theta)
  # Analytical: d/dmu = sum(x) - n*mu = 15 - 10 = 5
  expect_equal(s[1], sum_x - n * theta[1], tolerance = tol_mle)

  H <- hessian(loglik, theta)
  # Analytical: d^2/dmu^2 = -n
  expect_equal(H[1,1], -n, tolerance = tol_mle)

  I <- observed_information(loglik, theta)
  expect_equal(I[1,1], n, tolerance = tol_mle)
})

# -- Normal distribution (mu and sigma) ----------------------------------------

test_that("Normal (mu, sigma): score and hessian", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  # sum((x-mu)^2) = sum(x^2) - 2*mu*sum(x) + n*mu^2
  loglik <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2  # sum of squared residuals
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  mu0 <- 3
  sigma0 <- sqrt(2)
  theta <- c(mu0, sigma0)

  s <- score(loglik, theta)
  H <- hessian(loglik, theta)

  # Analytical score:
  ss0 <- sum((data - mu0)^2)
  expected_s1 <- sum(data - mu0) / sigma0^2
  expected_s2 <- -n / sigma0 + ss0 / sigma0^3
  expect_equal(s[1], expected_s1, tolerance = tol_mle)
  expect_equal(s[2], expected_s2, tolerance = tol_mle)

  # Analytical Hessian:
  expect_equal(H[1,1], -n/sigma0^2, tolerance = tol_mle)
  expect_equal(H[1,2], -2*sum(data-mu0)/sigma0^3, tolerance = tol_mle)
  expect_equal(H[2,1], H[1,2], tolerance = 1e-10)
  expect_equal(H[2,2], n/sigma0^2 - 3*ss0/sigma0^4, tolerance = tol_mle)
})

# -- Poisson distribution (lambda) --------------------------------------------

test_that("Poisson (lambda): score and hessian", {
  data <- c(1, 3, 2, 0, 4, 2)
  n <- length(data)
  sum_x <- sum(data)

  loglik <- function(theta) {
    lambda <- theta[1]
    sum_x * log(lambda) - n * lambda - sum(lfactorial(data))
  }

  lambda0 <- 2
  theta <- c(lambda0)

  s <- score(loglik, theta)
  H <- hessian(loglik, theta)

  expect_equal(s[1], sum_x/lambda0 - n, tolerance = tol_mle)
  expect_equal(H[1,1], -sum_x/lambda0^2, tolerance = tol_mle)
})

# -- Gamma distribution (shape, rate) -----------------------------------------

test_that("Gamma (shape, rate): score matches numerical", {
  data <- c(1.5, 2.1, 0.8, 3.2, 1.9)
  n <- length(data)
  sum_x <- sum(data)
  sum_log_x <- sum(log(data))

  loglik <- function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    n * alpha * log(beta) - n * lgamma(alpha) +
      (alpha - 1) * sum_log_x - beta * sum_x
  }

  theta <- c(2, 1.5)

  s <- score(loglik, theta)
  H <- hessian(loglik, theta)

  # Compare against numerical derivatives of the same function
  num_loglik <- function(t) {
    n * t[1] * log(t[2]) - n * lgamma(t[1]) +
      (t[1] - 1) * sum_log_x - t[2] * sum_x
  }
  num_s <- numerical_gradient(num_loglik, theta)
  num_H <- numerical_hessian(num_loglik, theta)

  expect_equal(s, num_s, tolerance = tol_mle)
  expect_equal(H, num_H, tolerance = 1e-4)
})

# -- score_and_hessian ---------------------------------------------------------

test_that("score_and_hessian from analytical score", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  # Analytical score for Normal(mu, sigma)
  score_fn <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    s1 <- (sum_x - n * mu) / sigma^2
    s2 <- -n / sigma + ss / sigma^3
    list(s1, s2)
  }

  mu0 <- 3
  sigma0 <- sqrt(2)
  theta <- c(mu0, sigma0)

  result <- score_and_hessian(score_fn, theta)

  expect_equal(result$score[1], sum(data - mu0) / sigma0^2, tolerance = tol_mle)
  expect_equal(result$score[2], -n/sigma0 + sum((data-mu0)^2)/sigma0^3,
               tolerance = tol_mle)

  expect_equal(result$hessian[1,1], -n/sigma0^2, tolerance = tol_mle)
  expect_equal(result$hessian[1,2], result$hessian[2,1], tolerance = 1e-10)
})

# -- Verify score_and_hessian Hessian matches hessian() from loglik -----------

test_that("score_and_hessian Hessian matches hessian() from loglik", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  loglik <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  score_fn <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    s1 <- (sum_x - n * mu) / sigma^2
    s2 <- -n / sigma + ss / sigma^3
    list(s1, s2)
  }

  theta <- c(3, sqrt(2))

  H_from_loglik <- hessian(loglik, theta)
  result <- score_and_hessian(score_fn, theta)

  expect_equal(result$hessian, H_from_loglik, tolerance = 1e-6)
})

# -- Works with optim-style function (vector indexing) -------------------------

test_that("MLE helpers work with optim-style functions", {
  data <- c(2.1, 3.5, 2.8, 4.0, 3.2)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  loglik <- function(theta) {
    mu <- theta[1]
    -0.5 * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  s <- score(loglik, c(3))
  expect_true(is.numeric(s))
  expect_equal(length(s), 1)

  H <- hessian(loglik, c(3))
  expect_true(is.matrix(H))
  expect_equal(dim(H), c(1, 1))
  expect_equal(H[1,1], -n, tolerance = tol_mle)
})

# -- Multi-parameter indexing --------------------------------------------------

test_that("multi-parameter loglik with vector indexing", {
  loglik <- function(theta) {
    x <- theta[1]
    y <- theta[2]
    -(x - 1)^2 - 2 * (y - 3)^2
  }

  theta <- c(0, 0)
  s <- score(loglik, theta)
  # d/dx = -2*(x-1) = 2 at x=0
  # d/dy = -4*(y-3) = 12 at y=0
  expect_equal(s[1], 2, tolerance = tol_mle)
  expect_equal(s[2], 12, tolerance = tol_mle)

  H <- hessian(loglik, theta)
  expect_equal(H[1,1], -2, tolerance = tol_mle)
  expect_equal(H[2,2], -4, tolerance = tol_mle)
  expect_equal(H[1,2], 0, tolerance = tol_mle)
  expect_equal(H[2,1], 0, tolerance = tol_mle)
})

# -- Exponential family --------------------------------------------------------

test_that("Exponential distribution: score and hessian", {
  data <- c(0.5, 1.2, 0.3, 2.1, 0.8)
  n <- length(data)
  sum_x <- sum(data)

  loglik <- function(theta) {
    rate <- theta[1]
    n * log(rate) - rate * sum_x
  }

  rate0 <- 1.5
  theta <- c(rate0)

  s <- score(loglik, theta)
  H <- hessian(loglik, theta)

  expect_equal(s[1], n/rate0 - sum_x, tolerance = tol_mle)
  expect_equal(H[1,1], -n/rate0^2, tolerance = tol_mle)
})
