# Tests for integration of score() and hessian() with R's optimizers.

# ============================================================================
# score() with optim()
# ============================================================================

test_that("score() works as optim() gr argument on a quadratic", {
  # Minimize f(x) = (x - 3)^2: minimum at x = 3
  f_ll <- function(theta) {
    x <- theta[1]
    -(x - 3)^2  # negative because we'll negate for minimization
  }

  nll <- function(theta) -f_ll(theta)
  ngr <- function(theta) -score(f_ll, theta)

  fit <- optim(0, fn = nll, gr = ngr, method = "BFGS")
  expect_equal(fit$par, 3, tolerance = 1e-6)
  expect_equal(fit$convergence, 0)
})

test_that("score() with optim() on 2D quadratic", {
  # Minimize f(x,y) = (x-1)^2 + (y-2)^2
  f_ll <- function(theta) {
    -(theta[1] - 1)^2 - (theta[2] - 2)^2
  }

  nll <- function(theta) -f_ll(theta)
  ngr <- function(theta) -score(f_ll, theta)

  fit <- optim(c(0, 0), fn = nll, gr = ngr, method = "BFGS")
  expect_equal(fit$par, c(1, 2), tolerance = 1e-6)
})

test_that("negated score gives correct minimization direction", {
  # score() gives the gradient of the log-likelihood (ascent direction)
  # -score() gives the descent direction for minimization
  f_ll <- function(theta) {
    mu <- theta[1]
    -0.5 * (5 - mu)^2
  }

  s <- score(f_ll, 3)  # d/dmu = (5 - 3) = 2 (positive: ascend toward 5)
  expect_equal(s[1], 2)

  neg_s <- -score(f_ll, 3)
  expect_equal(neg_s[1], -2)  # negative: descend toward 5 from below
})

# ============================================================================
# score() + hessian() with nlminb()
# ============================================================================

test_that("score() and hessian() work with nlminb()", {
  f_ll <- function(theta) {
    -(theta[1] - 3)^2 - 2 * (theta[2] + 1)^2
  }

  nll <- function(theta) -f_ll(theta)
  ngr <- function(theta) -score(f_ll, theta)
  nhess <- function(theta) -hessian(f_ll, theta)

  fit <- nlminb(c(0, 0), objective = nll, gradient = ngr, hessian = nhess)
  expect_equal(fit$par, c(3, -1), tolerance = 1e-6)
  expect_equal(fit$convergence, 0)
})

# ============================================================================
# AD gradient matches optim's numerical gradient (same solution)
# ============================================================================

test_that("optim with AD gradient matches optim without gradient", {
  set.seed(42)
  data_norm <- rnorm(50, mean = 5, sd = 2)
  n <- length(data_norm)
  sum_x <- sum(data_norm)
  sum_x2 <- sum(data_norm^2)

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  nll <- function(theta) -ll(theta)
  ngr <- function(theta) -score(ll, theta)

  # Use L-BFGS-B with lower bound on sigma to prevent log(negative) NaN warnings
  start <- c(3, 1.5)
  fit_no_gr <- optim(start, fn = nll, method = "L-BFGS-B",
                     lower = c(-Inf, 0.01))
  fit_ad_gr <- optim(start, fn = nll, gr = ngr, method = "L-BFGS-B",
                     lower = c(-Inf, 0.01))

  # Both should converge to the same point
  expect_equal(fit_ad_gr$par, fit_no_gr$par, tolerance = 1e-5)
  expect_equal(fit_ad_gr$convergence, 0)
})

# ============================================================================
# Standard errors from observed_information() match analytical values
# ============================================================================

test_that("SEs from observed_information match analytical for Normal(mu, sigma)", {
  fix <- make_normal_fixture()
  n <- fix$n; sum_x <- fix$sum_x; sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  # Find MLE
  mle <- c(fix$mle_mu, fix$mle_sigma)

  # Observed information at MLE
  obs_info <- observed_information(ll, mle)
  vcov <- solve(obs_info)
  se_ad <- sqrt(diag(vcov))

  # Analytical SEs at the MLE:
  # Var(mu_hat) = sigma^2/n, Var(sigma_hat) = sigma^2/(2n)
  se_analytical <- c(
    fix$mle_sigma / sqrt(n),
    fix$mle_sigma / sqrt(2 * n)
  )

  expect_equal(se_ad, se_analytical, tolerance = 1e-6)
})

# ============================================================================
# Poisson MLE via optim() matches sum(x)/n
# ============================================================================

test_that("Poisson MLE via optim + AD matches analytical", {
  set.seed(123)
  data_pois <- rpois(100, lambda = 3.5)
  n_pois <- length(data_pois)
  sum_x_pois <- sum(data_pois)
  sum_lfact <- sum(lfactorial(data_pois))

  ll <- function(theta) {
    lambda <- theta[1]
    sum_x_pois * log(lambda) - n_pois * lambda - sum_lfact
  }

  nll <- function(theta) -ll(theta)
  ngr <- function(theta) -score(ll, theta)

  fit <- optim(1, fn = nll, gr = ngr, method = "BFGS")

  analytical_mle <- sum_x_pois / n_pois
  expect_equal(fit$par, analytical_mle, tolerance = 1e-6)
})

# ============================================================================
# Logistic regression: optim + AD matches glm()
# ============================================================================

test_that("logistic regression via optim + AD matches glm()", {
  set.seed(7)
  n_lr <- 200
  x1 <- rnorm(n_lr)
  x2 <- rnorm(n_lr)
  X <- cbind(1, x1, x2)
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

  nll <- function(theta) -ll_num(theta)
  ngr <- function(theta) -score(ll_dual, theta)

  fit_optim <- optim(c(0, 0, 0), fn = nll, gr = ngr, method = "BFGS")
  fit_glm <- glm(y ~ x1 + x2, family = binomial)

  expect_equal(fit_optim$par, unname(coef(fit_glm)), tolerance = 1e-4)
})

test_that("logistic regression SEs from AD match glm()", {
  set.seed(7)
  n_lr <- 200
  x1 <- rnorm(n_lr)
  x2 <- rnorm(n_lr)
  X <- cbind(1, x1, x2)
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

  nll <- function(theta) -ll_num(theta)
  ngr <- function(theta) -score(ll_dual, theta)

  fit_optim <- optim(c(0, 0, 0), fn = nll, gr = ngr, method = "BFGS")
  fit_glm <- glm(y ~ x1 + x2, family = binomial)

  # SEs from AD Hessian
  obs_info <- observed_information(ll_dual, fit_optim$par)
  vcov_ad <- solve(obs_info)
  se_ad <- sqrt(diag(vcov_ad))

  se_glm <- summary(fit_glm)$coefficients[, "Std. Error"]

  expect_equal(se_ad, unname(se_glm), tolerance = 0.01)
})

# ============================================================================
# nlminb converges for multi-parameter model
# ============================================================================

test_that("nlminb with full AD converges for Normal(mu, sigma)", {
  fix <- make_normal_fixture()
  n <- fix$n; sum_x <- fix$sum_x; sum_x2 <- fix$sum_x2

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  nll <- function(theta) -ll(theta)
  ngr <- function(theta) -score(ll, theta)
  nhess <- function(theta) -hessian(ll, theta)

  fit <- nlminb(c(0, 1), objective = nll, gradient = ngr, hessian = nhess)

  expect_equal(fit$par[1], fix$mle_mu, tolerance = 1e-6)
  expect_equal(fit$par[2], fix$mle_sigma, tolerance = 1e-4)
  expect_equal(fit$convergence, 0)
})
