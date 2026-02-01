#' Arithmetic and comparison operators for dual numbers
#'
#' Implements arithmetic (\code{+}, \code{-}, \code{*}, \code{/}, \code{^}),
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
#' @aliases Ops,dualr,dualr-method Ops,dualr,numeric-method Ops,numeric,dualr-method
NULL

# -- Internal helper: recursive zero check for nested duals ------------------

.is_zero <- function(x) {
  if (identical(x, 0)) return(TRUE)
  if (is.numeric(x)) return(length(x) == 1L && x == 0)
  if (is(x, "dualr")) return(.is_zero(x@value) && .is_zero(x@deriv))
  FALSE
}

# =============================================================================
# Per-operation methods for the 5 hot-path arithmetic operators
# Each has 3 type combos: (dualr,dualr), (dualr,numeric), (numeric,dualr)
# =============================================================================

# -- Addition -----------------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("+", signature(e1 = "dualr", e2 = "dualr"), function(e1, e2) {
  .dual(e1@value + e2@value, e1@deriv + e2@deriv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("+", signature(e1 = "dualr", e2 = "numeric"), function(e1, e2) {
  .dual(e1@value + e2, e1@deriv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("+", signature(e1 = "numeric", e2 = "dualr"), function(e1, e2) {
  .dual(e1 + e2@value, e2@deriv)
})

# -- Subtraction --------------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("-", signature(e1 = "dualr", e2 = "dualr"), function(e1, e2) {
  .dual(e1@value - e2@value, e1@deriv - e2@deriv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("-", signature(e1 = "dualr", e2 = "numeric"), function(e1, e2) {
  .dual(e1@value - e2, e1@deriv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("-", signature(e1 = "numeric", e2 = "dualr"), function(e1, e2) {
  .dual(e1 - e2@value, -e2@deriv)
})

# -- Multiplication -----------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("*", signature(e1 = "dualr", e2 = "dualr"), function(e1, e2) {
  v1 <- e1@value; v2 <- e2@value
  .dual(v1 * v2, v1 * e2@deriv + v2 * e1@deriv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("*", signature(e1 = "dualr", e2 = "numeric"), function(e1, e2) {
  .dual(e1@value * e2, e1@deriv * e2)
})

#' @rdname dual-arithmetic
#' @export
setMethod("*", signature(e1 = "numeric", e2 = "dualr"), function(e1, e2) {
  .dual(e1 * e2@value, e1 * e2@deriv)
})

# -- Division -----------------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("/", signature(e1 = "dualr", e2 = "dualr"), function(e1, e2) {
  v2 <- e2@value
  .dual(e1@value / v2, (e1@deriv * v2 - e1@value * e2@deriv) / (v2 * v2))
})

#' @rdname dual-arithmetic
#' @export
setMethod("/", signature(e1 = "dualr", e2 = "numeric"), function(e1, e2) {
  .dual(e1@value / e2, e1@deriv / e2)
})

#' @rdname dual-arithmetic
#' @export
setMethod("/", signature(e1 = "numeric", e2 = "dualr"), function(e1, e2) {
  v2 <- e2@value
  .dual(e1 / v2, -(e1 * e2@deriv) / (v2 * v2))
})

# -- Power --------------------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("^", signature(e1 = "dualr", e2 = "dualr"), function(e1, e2) {
  v1 <- e1@value; d1 <- e1@deriv
  v2 <- e2@value; d2 <- e2@deriv
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
  .dual(val, drv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("^", signature(e1 = "dualr", e2 = "numeric"), function(e1, e2) {
  v1 <- e1@value
  .dual(v1^e2, e2 * v1^(e2 - 1) * e1@deriv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("^", signature(e1 = "numeric", e2 = "dualr"), function(e1, e2) {
  val <- e1^e2@value
  .dual(val, val * log(e1) * e2@deriv)
})

# =============================================================================
# Ops group generic fallback for remaining operators
# Handles: %%, %/%, ==, !=, <, >, <=, >=, &, |
# =============================================================================

#' @rdname dual-arithmetic
#' @export
setMethod("Ops", signature(e1 = "dualr", e2 = "dualr"), function(e1, e2) {
  v1 <- e1@value; d1 <- e1@deriv
  v2 <- e2@value; d2 <- e2@deriv

  switch(.Generic,
    "%%" = .dual(v1 %% v2, d1),
    "%/%" = .dual(v1 %/% v2, 0),

    "==" = v1 == v2,
    "!=" = v1 != v2,
    "<"  = v1 <  v2,
    ">"  = v1 >  v2,
    "<=" = v1 <= v2,
    ">=" = v1 >= v2,

    "&" = v1 & v2,
    "|" = v1 | v2,

    stop(sprintf("operator '%s' not implemented for dual numbers", .Generic))
  )
})

#' @rdname dual-arithmetic
#' @export
setMethod("Ops", signature(e1 = "dualr", e2 = "numeric"), function(e1, e2) {
  v1 <- e1@value

  switch(.Generic,
    "%%" = .dual(v1 %% e2, e1@deriv),
    "%/%" = .dual(v1 %/% e2, 0),

    "==" = v1 == e2,
    "!=" = v1 != e2,
    "<"  = v1 <  e2,
    ">"  = v1 >  e2,
    "<=" = v1 <= e2,
    ">=" = v1 >= e2,

    "&" = v1 & e2,
    "|" = v1 | e2,

    stop(sprintf("operator '%s' not implemented for dual numbers", .Generic))
  )
})

#' @rdname dual-arithmetic
#' @export
setMethod("Ops", signature(e1 = "numeric", e2 = "dualr"), function(e1, e2) {
  v2 <- e2@value

  switch(.Generic,
    "%%" = .dual(e1 %% v2, 0),
    "%/%" = .dual(e1 %/% v2, 0),

    "==" = e1 == v2,
    "!=" = e1 != v2,
    "<"  = e1 <  v2,
    ">"  = e1 >  v2,
    "<=" = e1 <= v2,
    ">=" = e1 >= v2,

    "&" = e1 & v2,
    "|" = e1 | v2,

    stop(sprintf("operator '%s' not implemented for dual numbers", .Generic))
  )
})

# -- Unary operators -----------------------------------------------------------

#' @rdname dual-arithmetic
#' @export
setMethod("+", signature(e1 = "dualr", e2 = "missing"), function(e1, e2) e1)

#' @rdname dual-arithmetic
#' @export
setMethod("-", signature(e1 = "dualr", e2 = "missing"), function(e1, e2) {
  .dual(-e1@value, -e1@deriv)
})

#' @rdname dual-arithmetic
#' @export
setMethod("!", signature(x = "dualr"), function(x) !x@value)

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
#' @aliases Summary,dualr-method
NULL

#' @rdname dual-summary
#' @export
setMethod("Summary", "dualr", function(x, ..., na.rm = FALSE) {
  args <- list(x, ...)

  switch(.Generic,
    "sum" = {
      args <- lapply(args, .as_dual)
      total_v <- 0
      total_d <- 0
      for (a in args) {
        total_v <- total_v + a@value
        total_d <- total_d + a@deriv
      }
      .dual(total_v, total_d)
    },
    "prod" = {
      args <- lapply(args, .as_dual)
      Reduce("*", args)
    },
    "min" = {
      args <- lapply(args, .as_dual)
      Reduce(.dual_min, args)
    },
    "max" = {
      args <- lapply(args, .as_dual)
      Reduce(.dual_max, args)
    },
    "range" = {
      args <- lapply(args, .as_dual)
      mn <- Reduce(.dual_min, args)
      mx <- Reduce(.dual_max, args)
      dual_vector(list(mn, mx))
    },
    "any" = any(vapply(args, function(a) as.logical(if (is(a, "dualr")) a@value else a), logical(1))),
    "all" = all(vapply(args, function(a) as.logical(if (is(a, "dualr")) a@value else a), logical(1))),
    stop(sprintf("Summary function '%s' not implemented for dual numbers", .Generic))
  )
})
