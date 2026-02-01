#' Math group generic for dual numbers
#'
#' @description
#' Implements all standard mathematical functions for dual numbers via
#' the chain rule: \code{f(dual(a, b)) = dual(f(a), df(a) * b)}.
#'
#' Supported functions: \code{abs}, \code{sign}, \code{sqrt}, \code{floor},
#' \code{ceiling}, \code{trunc}, \code{round},
#' \code{exp}, \code{expm1}, \code{log}, \code{log2}, \code{log10}, \code{log1p},
#' \code{cos}, \code{sin}, \code{tan}, \code{cospi}, \code{sinpi}, \code{tanpi},
#' \code{acos}, \code{asin}, \code{atan},
#' \code{cosh}, \code{sinh}, \code{tanh}, \code{acosh}, \code{asinh}, \code{atanh},
#' \code{gamma}, \code{lgamma}, \code{digamma}, \code{trigamma},
#' \code{cumsum}, \code{factorial}, \code{lfactorial}.
#'
#' @param x A \code{dual} number.
#' @return A \code{dual} with the function applied to the value and the
#'   derivative propagated via the chain rule.
#'
#' @examples
#' x <- dual_variable(pi / 4)
#' value(sin(x))  # sin(pi/4)
#' deriv(sin(x))  # cos(pi/4)
#'
#' y <- dual_variable(2)
#' value(exp(y))  # exp(2)
#' deriv(exp(y))  # exp(2)
#' deriv(log(y))  # 1/2
#'
#' @name dual-math
#' @aliases Math,dualr-method
NULL

# =============================================================================
# Per-function methods for the 3 hottest math functions: exp, log, sqrt
# These bypass the Math group generic + switch dispatch entirely.
# =============================================================================

#' @rdname dual-math
#' @export
setMethod("exp", "dualr", function(x) {
  ev <- exp(x@value)
  .dual(ev, x@deriv * ev)
})

#' @rdname dual-math
#' @export
setMethod("sqrt", "dualr", function(x) {
  sv <- sqrt(x@value)
  .dual(sv, x@deriv / (2 * sv))
})

# =============================================================================
# Math group generic fallback for remaining functions
# =============================================================================

#' @rdname dual-math
#' @export
setMethod("Math", "dualr", function(x) {
  v <- x@value
  d <- x@deriv

  switch(.Generic,
    # -- Absolute value / sign --
    "abs"     = .dual(abs(v), d * sign(v)),
    "sign"    = .dual(sign(v), 0),

    # -- Rounding (piecewise constant, derivative 0 almost everywhere) --
    "floor"   = .dual(floor(v), 0),
    "ceiling" = .dual(ceiling(v), 0),
    "trunc"   = .dual(trunc(v), 0),
    "round"   = .dual(round(v), 0),

    # Note: sqrt, exp, and log have dedicated setMethod() dispatches above
    # (and for log, below) that take priority over this group generic.

    # -- Exponential / logarithm --
    "expm1"   = .dual(expm1(v), d * exp(v)),
    "log2"    = .dual(log2(v), d / (v * log(2))),
    "log10"   = .dual(log10(v), d / (v * log(10))),
    "log1p"   = .dual(log1p(v), d / (1 + v)),

    # -- Trigonometric --
    "cos"     = .dual(cos(v), -d * sin(v)),
    "sin"     = .dual(sin(v), d * cos(v)),
    "tan"     = { cv <- cos(v); .dual(tan(v), d / (cv * cv)) },
    "cospi"   = .dual(cospi(v), -d * pi * sinpi(v)),
    "sinpi"   = .dual(sinpi(v), d * pi * cospi(v)),
    "tanpi"   = { cv <- cospi(v); .dual(tanpi(v), d * pi / (cv * cv)) },

    # -- Inverse trigonometric --
    "acos"    = .dual(acos(v), -d / sqrt(1 - v * v)),
    "asin"    = .dual(asin(v), d / sqrt(1 - v * v)),
    "atan"    = .dual(atan(v), d / (1 + v * v)),

    # -- Hyperbolic --
    "cosh"    = .dual(cosh(v), d * sinh(v)),
    "sinh"    = .dual(sinh(v), d * cosh(v)),
    "tanh"    = { tv <- tanh(v); .dual(tv, d * (1 - tv * tv)) },

    # -- Inverse hyperbolic --
    "acosh"   = .dual(acosh(v), d / sqrt(v * v - 1)),
    "asinh"   = .dual(asinh(v), d / sqrt(v * v + 1)),
    "atanh"   = .dual(atanh(v), d / (1 - v * v)),

    # -- Gamma-related (chain rule applied) --
    "gamma"   = {
      gv <- gamma(v)
      .dual(gv, d * gv * digamma(v))
    },
    "lgamma"  = .dual(lgamma(v), d * digamma(v)),
    "digamma" = .dual(digamma(v), d * trigamma(v)),
    "trigamma" = .dual(trigamma(v), d * psigamma(v, deriv = 2L)),

    # -- Cumulative / factorial (not differentiable in usual sense) --
    "cummax"  = .dual(cummax(v), 0),
    "cummin"  = .dual(cummin(v), 0),
    "cumsum"  = .dual(cumsum(v), cumsum(d)),
    "cumprod" = {
      stop("cumprod() is not supported for dual numbers (derivative cannot be propagated). Use Reduce(\"*\", ...) instead.")
    },
    "factorial" = {
      gv <- gamma(v + 1)
      .dual(gv, d * gv * digamma(v + 1))
    },
    "lfactorial" = .dual(lfactorial(v), d * digamma(v + 1)),

    stop(sprintf("Math function '%s' not implemented for dual numbers", .Generic))
  )
})

# -- Math2 group generic (round, signif) ----------------------------------------

#' Math2 group generic for dual numbers
#'
#' Implements \code{round} and \code{signif} for dual numbers. These are
#' piecewise constant functions, so the derivative is zero almost everywhere.
#'
#' @param x A \code{dual} number.
#' @param digits Integer; number of digits for rounding.
#' @return A \code{dual} with the rounded value and zero derivative.
#'
#' @examples
#' x <- dual_variable(3.14159)
#' value(round(x, 2))  # 3.14
#' deriv(round(x, 2))  # 0 (piecewise constant)
#'
#' @name dual-math2
#' @aliases Math2,dualr-method
NULL

#' @rdname dual-math2
#' @export
setMethod("Math2", "dualr", function(x, digits) {
  v <- x@value
  switch(.Generic,
    "round"  = .dual(round(v, digits), 0),
    "signif" = .dual(signif(v, digits), 0),
    stop(sprintf("Math2 function '%s' not implemented for dual numbers", .Generic))
  )
})

# -- Standalone math methods (not in Math group) --------------------------------

#' Two-argument arctangent for dual numbers
#'
#' @param y A dual or numeric.
#' @param x A dual or numeric.
#' @return A dual representing atan2(y, x) with correct derivative.
#' @examples
#' y <- dual_variable(1)
#' x <- dual_constant(1)
#' result <- atan2(y, x)
#' value(result)  # pi/4
#'
#' @name dual-atan2
#' @aliases atan2,dualr,dualr-method atan2,dualr,numeric-method atan2,numeric,dualr-method
NULL

#' @rdname dual-atan2
#' @export
setMethod("atan2", signature(y = "dualr", x = "dualr"), function(y, x) {
  vy <- y@value; dy <- y@deriv
  vx <- x@value; dx <- x@deriv
  denom <- vx * vx + vy * vy
  .dual(atan2(vy, vx), (vx * dy - vy * dx) / denom)
})

#' @rdname dual-atan2
#' @export
setMethod("atan2", signature(y = "dualr", x = "numeric"), function(y, x) {
  vy <- y@value; dy <- y@deriv
  denom <- x * x + vy * vy
  .dual(atan2(vy, x), (x * dy) / denom)
})

#' @rdname dual-atan2
#' @export
setMethod("atan2", signature(y = "numeric", x = "dualr"), function(y, x) {
  vx <- x@value; dx <- x@deriv
  denom <- vx * vx + y * y
  .dual(atan2(y, vx), (-y * dx) / denom)
})

# -- max and min ---------------------------------------------------------------

#' Piecewise max and min for dual numbers
#'
#' Compares on value and propagates the derivative of the selected branch.
#'
#' @param x A dual number.
#' @param ... Additional dual or numeric values.
#' @param na.rm Ignored.
#' @return A \code{dual} representing the max or min.
#' @examples
#' x <- dual_variable(3)
#' y <- dual_variable(5)
#' value(max(x, y))  # 5
#' value(min(x, y))  # 3
#'
#' @name dual-maxmin
#' @aliases max,dualr-method min,dualr-method
NULL

#' @rdname dual-maxmin
#' @export
setMethod("max", signature(x = "dualr"), function(x, ..., na.rm = FALSE) {
  args <- list(x, ...)
  if (length(args) == 1L) return(x)
  args <- lapply(args, .as_dual)
  Reduce(.dual_max, args)
})

#' @rdname dual-maxmin
#' @export
setMethod("min", signature(x = "dualr"), function(x, ..., na.rm = FALSE) {
  args <- list(x, ...)
  if (length(args) == 1L) return(x)
  args <- lapply(args, .as_dual)
  Reduce(.dual_min, args)
})

# -- log with base argument ---------------------------------------------------

#' Logarithm with optional base for dual numbers
#'
#' @param x A dual number.
#' @param base Numeric base (default: \code{exp(1)} for natural log).
#' @return A \code{dual} representing \code{log(x, base)}.
#' @examples
#' x <- dual_variable(8)
#' value(log(x, base = 2))  # 3
#' deriv(log(x, base = 2))  # 1 / (8 * log(2))
#'
#' @name dual-log
#' @aliases log,dualr-method
NULL

#' @rdname dual-log
#' @export
setMethod("log", signature(x = "dualr"), function(x, base = exp(1)) {
  v <- x@value
  d <- x@deriv
  if (missing(base) || identical(base, exp(1))) {
    .dual(log(v), d / v)
  } else {
    .dual(log(v, base), d / (v * log(base)))
  }
})
