# Finite difference utilities for verifying AD derivatives.
# These are test-only helpers, not exported.

#' Central difference first derivative estimate
#' @param f Function of one numeric argument
#' @param x Point at which to estimate f'(x)
#' @param h Step size (default 1e-7 for O(h^2) accuracy)
#' @return Numeric estimate of f'(x)
central_difference <- function(f, x, h = 1e-7) {
  (f(x + h) - f(x - h)) / (2 * h)
}

#' Numerical gradient via central differences
#' @param f Function of a numeric vector
#' @param x Point at which to estimate the gradient
#' @param h Step size
#' @return Numeric vector (gradient)
numerical_gradient <- function(f, x, h = 1e-7) {
  p <- length(x)
  grad <- numeric(p)
  for (i in seq_len(p)) {
    x_plus <- x_minus <- x
    x_plus[i] <- x[i] + h
    x_minus[i] <- x[i] - h
    grad[i] <- (f(x_plus) - f(x_minus)) / (2 * h)
  }
  grad
}

#' Numerical Hessian via central differences of gradient
#' @param f Function of a numeric vector
#' @param x Point at which to estimate the Hessian
#' @param h Step size (slightly larger than gradient's h for stability)
#' @return p x p numeric matrix
numerical_hessian <- function(f, x, h = 1e-5) {
  p <- length(x)
  H <- matrix(0, nrow = p, ncol = p)
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      x_pp <- x_pm <- x_mp <- x_mm <- x
      x_pp[i] <- x_pp[i] + h; x_pp[j] <- x_pp[j] + h
      x_pm[i] <- x_pm[i] + h; x_pm[j] <- x_pm[j] - h
      x_mp[i] <- x_mp[i] - h; x_mp[j] <- x_mp[j] + h
      x_mm[i] <- x_mm[i] - h; x_mm[j] <- x_mm[j] - h
      H[i, j] <- (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4 * h * h)
      H[j, i] <- H[i, j]
    }
  }
  H
}

#' Numerical second derivative via central differences
#' @param f Function of one numeric argument
#' @param x Point at which to estimate f''(x)
#' @param h Step size
#' @return Numeric estimate of f''(x)
numerical_second_deriv <- function(f, x, h = 1e-5) {
  (f(x + h) - 2 * f(x) + f(x - h)) / (h * h)
}

# -- Shared test fixture: Normal MLE data ------------------------------------

#' Generate a Normal MLE test fixture with seed 42
#' @return List with data, n, sum_x, sum_x2, xbar, mle_mu, mle_sigma
make_normal_fixture <- function() {
  set.seed(42)
  data_norm <- rnorm(100, mean = 5, sd = 2)
  n <- length(data_norm)
  sum_x <- sum(data_norm)
  sum_x2 <- sum(data_norm^2)
  xbar <- mean(data_norm)
  list(data = data_norm, n = n, sum_x = sum_x, sum_x2 = sum_x2,
       xbar = xbar, mle_mu = xbar,
       mle_sigma = sqrt(mean((data_norm - xbar)^2)))
}

# -- Newton-Raphson optimizer for tests --------------------------------------

#' Newton-Raphson optimizer using AD score and Hessian
#' @param loglik Log-likelihood function of theta vector
#' @param theta0 Initial parameter values
#' @param tol Convergence tolerance on max absolute score
#' @param max_iter Maximum iterations
#' @param trace If TRUE, record all iterates
#' @return List with estimate, iterations, and optionally trace matrix
newton_raphson <- function(loglik, theta0, tol = 1e-8, max_iter = 50,
                           trace = FALSE) {
  theta <- theta0
  trace_list <- if (trace) list(theta) else NULL
  for (iter in seq_len(max_iter)) {
    s <- score(loglik, theta)
    H <- hessian(loglik, theta)
    step <- solve(H, s)
    theta <- theta - step
    if (trace) trace_list[[iter + 1L]] <- theta
    if (max(abs(s)) < tol) break
  }
  result <- list(estimate = theta, iterations = iter)
  if (trace) result$trace <- do.call(rbind, trace_list)
  result
}
