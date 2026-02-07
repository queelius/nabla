# Multi-parameter derivative functions (gradient, Hessian, Jacobian, D)
#
# Convention: f(x) where x supports x[1], x[2], ... indexing.
# This works with both numeric vectors and dual_vector objects.

# -- Internal helpers ----------------------------------------------------------

.make_dual_vector <- function(x, seed_index) {
  p <- length(x)
  duals <- vector("list", p)
  for (j in seq_len(p)) {
    duals[[j]] <- if (j == seed_index) .dual(x[j], 1) else .dual(x[j], 0)
  }
  new("dual_vector", duals)
}

#' Compute the gradient of a scalar-valued function
#'
#' Evaluates the gradient of \code{f} at \code{x} using forward-mode AD.
#' Equivalent to \code{D(f, x)} for scalar-valued \code{f}.
#'
#' The function should use \code{x[1]}, \code{x[2]}, etc. to access
#' parameters. This is compatible with the standard R convention used by
#' \code{optim}.
#'
#' @param f A function taking a parameter vector (numeric or dual) and
#'   returning a scalar. Parameters are accessed via \code{x[i]}.
#' @param x A numeric vector of parameter values.
#' @return A numeric vector of length \code{p} containing the gradient.
#' @export
#' @examples
#' f <- function(x) -(x[1] - 3)^2 - (x[2] - 5)^2
#' gradient(f, c(1, 2))
gradient <- function(f, x) {
  D(f, x, order = 1L)
}

#' Compute the Hessian of a scalar-valued function
#'
#' Computes the matrix of second partial derivatives of \code{f} at \code{x}
#' using forward-mode AD. Equivalent to \code{D(f, x, order = 2)}.
#'
#' @param f A function taking a parameter vector and returning a scalar.
#' @param x A numeric vector of parameter values.
#' @return A \code{p x p} numeric matrix (the Hessian).
#' @export
#' @examples
#' f <- function(x) -(x[1] - 3)^2 - (x[2] - 5)^2
#' hessian(f, c(1, 2))
hessian <- function(f, x) {
  D(f, x, order = 2L)
}

#' Compute the Jacobian of a vector-valued function
#'
#' Computes the first derivative of \code{f} at \code{x} using forward-mode
#' AD. Equivalent to \code{D(f, x)}. For \code{f: R^p -> R^m}, returns an
#' \code{m x p} matrix. For scalar-valued \code{f}, returns a length-p
#' gradient vector.
#'
#' @param f A function taking a parameter vector and returning a scalar,
#'   list of scalars, or a \code{dual_vector}.
#' @param x A numeric vector of parameter values (length \code{p}).
#' @return An \code{m x p} numeric matrix for vector-valued \code{f}, or
#'   a numeric vector of length \code{p} for scalar-valued \code{f}.
#' @export
#' @examples
#' f <- function(x) {
#'   a <- x[1]; b <- x[2]
#'   list(a * b, a^2, sin(b))
#' }
#' jacobian(f, c(2, pi/4))
#'
#' g <- function(x) x[1]^2 + x[2]^2
#' jacobian(g, c(3, 4))
jacobian <- function(f, x) {
  D(f, x, order = 1L)
}

# == Composable total derivative operator D ====================================

# -- Internal tensor helpers ---------------------------------------------------

.D_output_shape <- function(result) {
  if (is(result, "dualr")) return(integer(0))             # scalar
  if (is.numeric(result) && length(result) == 1L) return(integer(0))
  if (is(result, "dual_vector")) return(length(result))
  if (is.list(result)) {
    d <- dim(result)
    if (!is.null(d)) return(d)
    return(length(result))
  }
  if (is.array(result)) return(dim(result))
  if (is.numeric(result)) return(length(result))
  stop("D: unsupported return type from f")
}

.D_flatten <- function(result) {
  if (is(result, "dualr")) return(list(result))
  if (is.numeric(result) && length(result) == 1L) return(list(result))
  if (is(result, "dual_vector")) return(result@.Data)
  if (is.list(result)) return(result)    # includes dim-attributed lists
  if (is.array(result) || is.numeric(result)) return(as.list(as.vector(result)))
  stop("D: unsupported return type from f")
}

.D_extract_deriv <- function(x) {
  if (is(x, "dualr")) return(x@deriv)
  0
}

# all_derivs: length m*n, column-major [j=1 block, j=2 block, ...]
.D_build_tensor <- function(all_derivs, out_shape, n) {
  tensor_shape <- if (length(out_shape) == 0L) n else c(out_shape, n)
  has_dual <- FALSE
  for (d in all_derivs) {
    if (is(d, "dualr")) { has_dual <- TRUE; break }
  }
  if (!has_dual) {
    data <- vapply(all_derivs, function(x) {
      if (is.numeric(x)) as.double(x) else 0
    }, numeric(1))
    if (length(tensor_shape) == 1L) return(data)
    return(array(data, dim = tensor_shape))
  }
  # List-array for dual elements (consumed by outer D)
  dim(all_derivs) <- tensor_shape
  all_derivs
}

# -- D operator ----------------------------------------------------------------

#' Composable total derivative operator
#'
#' \code{D(f)} returns the derivative of \code{f} as a new function.
#' \code{D(f, x)} evaluates the derivative at \code{x}.
#' \code{D(D(f))} composes for second-order derivatives, and so on.
#'
#' Each application of \code{D} appends one \code{n}-dimension to the output
#' shape, where \code{n = length(x)}:
#' \itemize{
#'   \item For \code{f: R^n -> R}: D gives \code{(n)} gradient, D^2 gives
#'     \code{(n,n)} Hessian, D^3 gives \code{(n,n,n)}, etc.
#'   \item For \code{f: R^n -> R^m}: D gives \code{(m,n)} Jacobian, D^2 gives
#'     \code{(m,n,n)}, etc.
#' }
#'
#' The composability works because the S4 dispatch for \code{dualr} arithmetic
#' handles nested duals recursively. When \code{D(f)} is called with dual
#' inputs (from an outer \code{D}), derivative propagation is automatic.
#'
#' \code{gradient()}, \code{hessian()}, and \code{jacobian()} are convenience
#' wrappers: \code{gradient(f, x)} is \code{D(f, x)}, \code{hessian(f, x)} is
#' \code{D(f, x, order = 2)}, and \code{jacobian(f, x)} is \code{D(f, x)}.
#'
#' @param f A function taking a parameter vector (via \code{x[i]} indexing)
#'   and returning a scalar, list, or \code{dual_vector}.
#' @param x Optional numeric vector. If provided, evaluates \code{D(f)(x)}.
#' @param order Derivative order (default 1). \code{order = k} applies
#'   \code{D} k times.
#' @return If \code{x} is \code{NULL}, a function. Otherwise, a numeric
#'   vector, matrix, or array of the appropriate tensor shape.
#' @export
#' @examples
#' f <- function(x) x[1]^2 * x[2]
#' D(f, c(3, 4))
#' D(f, c(3, 4), order = 2)
#'
#' g <- function(x) list(x[1] * x[2], x[1]^2)
#' D(g, c(2, 3))
#'
#' Df <- D(f)
#' DDf <- D(Df)
#' DDf(c(3, 4))
D <- function(f, x = NULL, order = 1L) {
  order <- as.integer(order)
  if (order < 1L) stop("order must be a positive integer")
  g <- f
  for (i in seq_len(order)) g <- .make_D_once(g)
  if (is.null(x)) return(g)
  g(x)
}

.make_D_once <- function(f) {
  force(f)
  function(x) .eval_D(f, x)
}

.eval_D <- function(f, x) {
  n <- length(x)
  x1 <- .make_dual_vector(x, 1L)
  result1 <- f(x1)
  out_shape <- .D_output_shape(result1)
  elems1 <- .D_flatten(result1)
  m <- length(elems1)
  all_derivs <- vector("list", m * n)
  for (k in seq_len(m)) all_derivs[[k]] <- .D_extract_deriv(elems1[[k]])
  if (n > 1L) {
    for (j in 2:n) {
      xj <- .make_dual_vector(x, j)
      elems_j <- .D_flatten(f(xj))
      offset <- (j - 1L) * m
      for (k in seq_len(m)) {
        all_derivs[[offset + k]] <- .D_extract_deriv(elems_j[[k]])
      }
    }
  }
  .D_build_tensor(all_derivs, out_shape, n)
}
