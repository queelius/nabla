#' @title dualr: Forward-Mode Automatic Differentiation via Dual Numbers
#'
#' @description
#' Implements forward-mode automatic differentiation using dual numbers
#' with S4 classes. Supports first- and second-order derivatives through
#' nested duals, with convenience functions for maximum likelihood
#' estimation workflows including score vectors, Hessian matrices, and
#' observed information.
#'
#' @section Core Types:
#' \describe{
#'   \item{\code{\link{dual}}}{Constructor for dual numbers.}
#'   \item{\code{\link{dual_variable}}}{Shorthand for \code{dual(x, 1)}.}
#'   \item{\code{\link{dual_constant}}}{Shorthand for \code{dual(x, 0)}.}
#'   \item{\code{\link{dual_vector}}}{Container for indexable dual vectors.}
#' }
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{\link{value}}}{Extract the primal value.}
#'   \item{\code{\link{deriv}}}{Extract the derivative component.}
#' }
#'
#' @section Second-Order:
#' \describe{
#'   \item{\code{\link{dual2_variable}}}{Nested dual for second derivatives.}
#'   \item{\code{\link{differentiate2}}}{Compute f(x), f'(x), f''(x) in one call.}
#' }
#'
#' @section MLE Helpers:
#' \describe{
#'   \item{\code{\link{score}}}{Gradient of a log-likelihood.}
#'   \item{\code{\link{hessian}}}{Hessian matrix of a log-likelihood.}
#'   \item{\code{\link{observed_information}}}{Negative Hessian.}
#'   \item{\code{\link{score_and_hessian}}}{Gradient + Hessian from a score function.}
#' }
#'
#' @references
#' Baydin, A. G., Pearlmutter, B. A., Radul, A. A., & Siskind, J. M. (2018).
#' Automatic differentiation in machine learning: a survey.
#' \emph{Journal of Machine Learning Research}, 18(153), 1--43.
#'
#' @seealso
#' Related CRAN packages: \pkg{dual}, \pkg{numDeriv}, \pkg{madness}
#'
#' @import methods
#' @importFrom stats pnorm
#' @docType package
#' @name dualr-package
#' @aliases dualr
"_PACKAGE"
