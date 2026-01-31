#' Arithmetic and comparison operators for dual numbers
#'
#' Implements the \code{Ops} group generic for dual numbers, supporting
#' arithmetic (\code{+}, \code{-}, \code{*}, \code{/}, \code{^}),
#' comparison (\code{==}, \code{!=}, \code{<}, \code{>}, \code{<=}, \code{>=}),
#' and logical (\code{&}, \code{|}) operators. Derivatives follow standard
#' calculus rules (sum, product, quotient, power, chain).
#'
#' @param e1,e2 Dual or numeric operands.
#' @param x A dual number (for unary \code{!}).
#' @return A \code{dual} for arithmetic ops; logical for comparisons.
#'
#' @examples
#' x <- dual_variable(3)
#' y <- dual_variable(4)
#'
#' # Arithmetic
#' value(x + y)   # 7
#' deriv(x * x)   # 6 (= 2 * 3)
#' value(x^2)     # 9
#' deriv(x^2)     # 6
#'
#' # Comparison (uses values only)
#' x < y   # TRUE
#' x == y  # FALSE
#'
#' @name dual-arithmetic
#' @aliases Ops,dual,dual-method Ops,dual,numeric-method Ops,numeric,dual-method
NULL

# -- Internal helper: recursive zero check for nested duals ------------------

.is_zero <- function(x) {
  if (identical(x, 0)) return(TRUE)
  if (is.numeric(x) && length(x) == 1L && x == 0) return(TRUE)
  if (is(x, "dual")) return(.is_zero(value(x)) && .is_zero(deriv(x)))
  FALSE
}

# -- dual, dual ---------------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("Ops", signature(e1 = "dual", e2 = "dual"), function(e1, e2) {
  v1 <- value(e1); d1 <- deriv(e1)
  v2 <- value(e2); d2 <- deriv(e2)

  switch(.Generic,
    "+" = dual(v1 + v2, d1 + d2),
    "-" = dual(v1 - v2, d1 - d2),
    "*" = dual(v1 * v2, v1 * d2 + v2 * d1),
    "/" = dual(v1 / v2, (d1 * v2 - v1 * d2) / (v2 * v2)),
    "^" = {
      val <- v1^v2
      is_d1_zero <- .is_zero(d1)
      is_d2_zero <- .is_zero(d2)
      if (is_d2_zero) {
        drv <- v2 * v1^(v2 - 1) * d1
      } else if (is_d1_zero) {
        drv <- val * log(v1) * d2
      } else {
        drv <- val * (d2 * log(v1) + v2 * d1 / v1)
      }
      dual(val, drv)
    },
    "%%" = dual(v1 %% v2, d1),
    "%/%" = dual(v1 %/% v2, 0),

    "==" = value(e1) == value(e2),
    "!=" = value(e1) != value(e2),
    "<"  = value(e1) <  value(e2),
    ">"  = value(e1) >  value(e2),
    "<=" = value(e1) <= value(e2),
    ">=" = value(e1) >= value(e2),

    "&" = value(e1) & value(e2),
    "|" = value(e1) | value(e2),

    stop(sprintf("operator '%s' not implemented for dual numbers", .Generic))
  )
})

# -- dual, numeric -------------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("Ops", signature(e1 = "dual", e2 = "numeric"), function(e1, e2) {
  callGeneric(e1, dual(e2, 0))
})

# -- numeric, dual -------------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("Ops", signature(e1 = "numeric", e2 = "dual"), function(e1, e2) {
  callGeneric(dual(e1, 0), e2)
})

# -- Unary operators -----------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("+", signature(e1 = "dual", e2 = "missing"), function(e1, e2) e1)

#' @rdname dual-arithmetic
#' @export
setMethod("-", signature(e1 = "dual", e2 = "missing"), function(e1, e2) {
  dual(-value(e1), -deriv(e1))
})

#' @rdname dual-arithmetic
#' @export
setMethod("!", signature(x = "dual"), function(x) !value(x))

# -- Summary group generic (sum, prod, min, max, range, any, all) -------------

#' Summary group generic for dual numbers
#'
#' Implements \code{sum}, \code{prod}, \code{min}, \code{max}, \code{range},
#' \code{any}, and \code{all} for dual numbers. Derivatives are propagated
#' correctly through \code{sum} (additive) and \code{prod} (multiplicative).
#'
#' @param x A dual number.
#' @param ... Additional dual or numeric values.
#' @param na.rm Logical; ignored (present for generic compatibility).
#' @return A \code{dual} for sum/prod/min/max; a \code{dual_vector} for range;
#'   logical for any/all.
#'
#' @examples
#' x <- dual_variable(2)
#' y <- dual_variable(5)
#' value(sum(x, y))   # 7
#' value(prod(x, y))  # 10
#'
#' @name dual-summary
#' @aliases Summary,dual-method
NULL

#' @rdname dual-summary
#' @export
setMethod("Summary", "dual", function(x, ..., na.rm = FALSE) {
  args <- list(x, ...)
  args <- lapply(args, .as_dual)

  switch(.Generic,
    "sum" = Reduce("+", args),
    "prod" = Reduce("*", args),
    "min" = Reduce(function(a, b) if (value(a) <= value(b)) a else b, args),
    "max" = Reduce(function(a, b) if (value(a) >= value(b)) a else b, args),
    "range" = {
      mn <- Reduce(function(a, b) if (value(a) <= value(b)) a else b, args)
      mx <- Reduce(function(a, b) if (value(a) >= value(b)) a else b, args)
      dual_vector(list(mn, mx))
    },
    "any" = any(sapply(args, function(a) as.logical(value(a)))),
    "all" = all(sapply(args, function(a) as.logical(value(a)))),
    stop(sprintf("Summary function '%s' not implemented for dual numbers", .Generic))
  )
})
