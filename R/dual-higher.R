# Higher-order derivatives via nested dual numbers
#
# A second-order dual encodes f(x), f'(x), and f''(x) simultaneously.
# The structure is: dual(dual(x, 1), dual(1, 0))
#
# After evaluating f on this structure:
#   value(value(result)) = f(x)
#   deriv(value(result)) = f'(x)
#   deriv(deriv(result))  = f''(x)
#
# The outer dual tracks d/dx, and the inner dual tracks the same variable.
# The cross-term in the outer deriv's deriv slot accumulates the second
# derivative via the product/chain rules.

#' Create a second-order dual variable
#'
#' Constructs a nested dual suitable for computing both first and second
#' derivatives of a function at a point.
#'
#' @param x A numeric value.
#' @return A nested dual: \code{dual(dual(x, 1), dual(1, 0))}.
#' @export
#' @examples
#' x <- dual2_variable(2)
#' result <- x^3
#' second_deriv(result)  # 12 = 6*x at x=2
dual2_variable <- function(x) {
  .dual(.dual(x, 1), .dual(1, 0))
}

#' Create a second-order dual constant
#'
#' Wraps a numeric as a nested dual with all derivative components zero.
#'
#' @param x A numeric value.
#' @return A nested dual: \code{dual(dual(x, 0), dual(0, 0))}.
#' @examples
#' k <- dual2_constant(5)
#' second_deriv(k^2)  # 0 (constant has no derivatives)
#' @export
dual2_constant <- function(x) {
  .dual(.dual(x, 0), .dual(0, 0))
}

# -- Extraction helpers --------------------------------------------------------

#' Extract function value from a second-order dual
#'
#' @param d A nested dual (result of evaluating a function on
#'   \code{dual2_variable}).
#' @return The numeric function value \eqn{f(x)}.
#' @examples
#' x <- dual2_variable(3)
#' result <- x^2
#' value2(result)  # 9
#' @export
value2 <- function(d) {
  d@value@value
}

#' Extract first derivative from a second-order dual
#'
#' @param d A nested dual.
#' @return The numeric first derivative \eqn{f'(x)}.
#' @examples
#' x <- dual2_variable(3)
#' result <- x^2
#' first_deriv(result)  # 6 (= 2*x at x=3)
#' @export
first_deriv <- function(d) {
  d@value@deriv
}

#' Extract second derivative from a second-order dual
#'
#' @param d A nested dual.
#' @return The numeric second derivative \eqn{f''(x)}.
#' @examples
#' x <- dual2_variable(3)
#' result <- x^2
#' second_deriv(result)  # 2
#' @export
second_deriv <- function(d) {
  d@deriv@deriv
}

# -- Convenience evaluator ----------------------------------------------------

#' Compute value, first and second derivatives of a function
#'
#' Evaluates \code{f} at a second-order dual variable seeded at \code{x},
#' returning all three quantities.
#'
#' @param f A function of one argument (numeric or dual).
#' @param x A numeric value at which to differentiate.
#' @return A named list with components \code{value}, \code{first}, and
#'   \code{second}.
#' @export
#' @examples
#' differentiate2(sin, pi/4)
#' # $value = sin(pi/4)
#' # $first = cos(pi/4)
#' # $second = -sin(pi/4)
differentiate2 <- function(f, x) {
  result <- f(dual2_variable(x))
  list(
    value  = result@value@value,
    first  = result@value@deriv,
    second = result@deriv@deriv
  )
}
