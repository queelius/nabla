# Special mathematical functions for dual numbers
#
# Functions not covered by the Math or Ops group generics.

# -- Error functions -----------------------------------------------------------

#' Error function
#'
#' Computes \eqn{\mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt}.
#'
#' @param x A numeric or dual value.
#' @return The error function value. For dual input, returns a dual with
#'   derivative \eqn{(2/\sqrt{\pi}) e^{-x^2}}.
#' @examples
#' erf(1)  # numeric
#' x <- dual_variable(1)
#' value(erf(x))  # erf(1)
#' @export
setGeneric("erf", function(x) standardGeneric("erf"))

#' @rdname erf
#' @export
setMethod("erf", "numeric", function(x) {
  2 * pnorm(x * sqrt(2)) - 1
})

#' @rdname erf
#' @export
setMethod("erf", "dual", function(x) {
  v <- value(x)
  d <- deriv(x)
  dual(erf(v), d * (2 / sqrt(pi)) * exp(-v * v))
})

#' Complementary error function
#'
#' Computes \eqn{\mathrm{erfc}(x) = 1 - \mathrm{erf}(x)}.
#'
#' @param x A numeric or dual value.
#' @return The complementary error function value.
#' @examples
#' erfc(1)  # 1 - erf(1)
#' x <- dual_variable(0)
#' value(erfc(x))  # 1
#' @export
setGeneric("erfc", function(x) standardGeneric("erfc"))

#' @rdname erfc
#' @export
setMethod("erfc", "numeric", function(x) {
  2 * pnorm(-x * sqrt(2))
})

#' @rdname erfc
#' @export
setMethod("erfc", "dual", function(x) {
  v <- value(x)
  d <- deriv(x)
  dual(erfc(v), -d * (2 / sqrt(pi)) * exp(-v * v))
})

# -- Beta functions ------------------------------------------------------------

#' Log-beta function for dual numbers
#'
#' @param a A numeric or dual value.
#' @param b A numeric or dual value.
#' @return \code{lbeta(a, b)} with derivative via digamma.
#' @examples
#' lbeta(2, 3)  # log(beta(2, 3))
#' a <- dual_variable(2)
#' value(lbeta(a, 3))  # log(beta(2, 3))
#' @export
setGeneric("lbeta", function(a, b) standardGeneric("lbeta"))

#' @rdname lbeta
#' @export
setMethod("lbeta", signature(a = "numeric", b = "numeric"), function(a, b) {
  base::lbeta(a, b)
})

#' @rdname lbeta
#' @export
setMethod("lbeta", signature(a = "dual", b = "dual"), function(a, b) {
  va <- value(a); da <- deriv(a)
  vb <- value(b); db <- deriv(b)
  lbv <- base::lbeta(va, vb)
  # d/da lbeta(a,b) = digamma(a) - digamma(a+b)
  # d/db lbeta(a,b) = digamma(b) - digamma(a+b)
  dab <- digamma(va + vb)
  drv <- da * (digamma(va) - dab) + db * (digamma(vb) - dab)
  dual(lbv, drv)
})

#' @rdname lbeta
#' @export
setMethod("lbeta", signature(a = "dual", b = "numeric"), function(a, b) {
  lbeta(a, dual(b, 0))
})

#' @rdname lbeta
#' @export
setMethod("lbeta", signature(a = "numeric", b = "dual"), function(a, b) {
  lbeta(dual(a, 0), b)
})

#' Beta function for dual numbers
#'
#' Computed as \code{exp(lbeta(a, b))}.
#'
#' @param a A numeric or dual value.
#' @param b A numeric or dual value.
#' @return \code{beta(a, b)} with derivative.
#' @examples
#' beta(2, 3)  # 1/12
#' a <- dual_variable(2)
#' value(beta(a, 3))
#' @export
setGeneric("beta", function(a, b) standardGeneric("beta"))

#' @rdname beta
#' @export
setMethod("beta", signature(a = "numeric", b = "numeric"), function(a, b) {
  base::beta(a, b)
})

#' @rdname beta
#' @export
setMethod("beta", signature(a = "dual", b = "dual"), function(a, b) {
  exp(lbeta(a, b))
})

#' @rdname beta
#' @export
setMethod("beta", signature(a = "dual", b = "numeric"), function(a, b) {
  exp(lbeta(a, dual(b, 0)))
})

#' @rdname beta
#' @export
setMethod("beta", signature(a = "numeric", b = "dual"), function(a, b) {
  exp(lbeta(dual(a, 0), b))
})

# -- Polygamma function --------------------------------------------------------

#' Polygamma function for dual numbers
#'
#' Computes \eqn{\psi^{(n)}(x)} where \eqn{n} is the derivative order.
#' The derivative of \eqn{\psi^{(n)}(x)} is \eqn{\psi^{(n+1)}(x)}.
#'
#' @param x A numeric or dual value.
#' @param deriv Integer derivative order (0 = digamma, 1 = trigamma, etc.).
#' @return The polygamma function value.
#' @examples
#' psigamma(1, deriv = 0)  # digamma(1)
#' x <- dual_variable(2)
#' value(psigamma(x, deriv = 1))  # trigamma(2)
#' @export
setGeneric("psigamma", function(x, deriv = 0L) standardGeneric("psigamma"))

#' @rdname psigamma
#' @export
setMethod("psigamma", signature(x = "numeric"), function(x, deriv = 0L) {
  base::psigamma(x, deriv = deriv)
})

#' @rdname psigamma
#' @export
setMethod("psigamma", signature(x = "dual"), function(x, deriv = 0L) {
  v <- value(x)
  d <- deriv(x)
  dual(base::psigamma(v, deriv = deriv),
       d * base::psigamma(v, deriv = deriv + 1L))
})
