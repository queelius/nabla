# Higher-order derivatives via nested dual numbers
#
# Seeding structure for order n:
#   n=0: just x (plain numeric)
#   n=1: dual(x, 1)
#   n=2: dual(dual(x, 1), dual(1, 0))
#   n=3: dual(dual(dual(x, 1), dual(1, 0)), dual(dual(1, 0), dual(0, 0)))

# =============================================================================
# Arbitrary-order constructors
# =============================================================================

#' Create a dual seeded for n-th order differentiation
#'
#' Recursively nests dual numbers to enable exact computation of
#' derivatives up to order \code{n}. The variable is seeded so that
#' after evaluating a function \code{f}, the k-th derivative can be
#' extracted with \code{deriv_n(result, k)}.
#'
#' @param x A numeric value at which to differentiate.
#' @param order A positive integer specifying the derivative order.
#' @return A (possibly nested) dual number.
#' @export
#' @examples
#' # 3rd derivative of x^4 at x=2: 4*3*2*x = 24*2 = 48
#' x <- dual_variable_n(2, order = 3)
#' r <- x^4
#' deriv_n(r, 3)  # 48
dual_variable_n <- function(x, order) {
  order <- as.integer(order)
  if (order < 0L) stop("order must be a non-negative integer")
  if (order == 0L) return(x)
  if (order == 1L) return(.dual(x, 1))
  .dual(dual_variable_n(x, order - 1L),
        dual_constant_n(1, order - 1L))
}

#' Create a constant dual for n-th order differentiation
#'
#' Wraps a numeric value as a nested dual with all derivative components
#' zero, representing a constant with respect to the differentiation
#' variable at nesting depth \code{n}.
#'
#' @param x A numeric value.
#' @param order A non-negative integer specifying the nesting depth.
#' @return A (possibly nested) dual number with zero derivatives.
#' @export
#' @examples
#' k <- dual_constant_n(5, order = 3)
#' deriv_n(k, 1)  # 0
#' deriv_n(k, 2)  # 0
#' deriv_n(k, 3)  # 0
dual_constant_n <- function(x, order) {
  order <- as.integer(order)
  if (order < 0L) stop("order must be a non-negative integer")
  if (order == 0L) return(x)
  if (order == 1L) return(.dual(x, 0))
  .dual(dual_constant_n(x, order - 1L),
        dual_constant_n(0, order - 1L))
}

# =============================================================================
# Extraction
# =============================================================================

#' Extract the k-th derivative from a nested dual result
#'
#' After evaluating a function on a dual created by
#' \code{\link{dual_variable_n}}, use \code{deriv_n} to extract any
#' derivative from 0 (the function value) up to the seeded order.
#'
#' @param d A (possibly nested) dual number, or a numeric.
#' @param k A non-negative integer: 0 for the function value, 1 for the
#'   first derivative, etc.
#' @return A numeric value.
#' @export
#' @examples
#' x <- dual_variable_n(1, order = 3)
#' r <- exp(x)
#' deriv_n(r, 0)  # exp(1) = 2.718...
#' deriv_n(r, 1)  # exp(1)
#' deriv_n(r, 2)  # exp(1)
#' deriv_n(r, 3)  # exp(1)
deriv_n <- function(d, k) {
  k <- as.integer(k)
  if (k < 0L) stop("k must be a non-negative integer")
  if (k == 0L) {
    while (is(d, "dualr")) d <- d@value
    return(d)
  }
  if (!is(d, "dualr")) return(0)
  deriv_n(d@deriv, k - 1L)
}

# =============================================================================
# Convenience evaluator
# =============================================================================

#' Compute a function value and all derivatives up to order n
#'
#' Evaluates \code{f} at a dual variable seeded for order \code{n},
#' returning the function value and all derivatives from 1 to \code{n}.
#'
#' @param f A function of one numeric argument.
#' @param x A numeric value at which to differentiate.
#' @param order A positive integer: the maximum derivative order.
#' @return A named list with components \code{value}, \code{d1},
#'   \code{d2}, ..., \code{d<order>}.
#' @export
#' @examples
#' # All derivatives of sin(x) at x = pi/4
#' differentiate_n(sin, pi/4, order = 4)
#' # $value = sin(pi/4)
#' # $d1 = cos(pi/4)
#' # $d2 = -sin(pi/4)
#' # $d3 = -cos(pi/4)
#' # $d4 = sin(pi/4)
differentiate_n <- function(f, x, order) {
  order <- as.integer(order)
  if (order < 1L) stop("order must be a positive integer")
  result <- f(dual_variable_n(x, order))
  out <- list(value = deriv_n(result, 0L))
  for (k in seq_len(order)) {
    out[[paste0("d", k)]] <- deriv_n(result, k)
  }
  out
}
