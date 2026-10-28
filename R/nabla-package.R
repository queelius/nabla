#' @title nabla: Exact Derivatives via Automatic Differentiation
#'
#' @description
#' Implements forward-mode automatic differentiation using dual numbers
#' with S4 classes. Supports exact arbitrary-order derivatives through
#' recursive nesting of duals, with high-level functions for computing
#' gradients, Hessian matrices, and Jacobians of functions composed of
#' supported scalar operations. Loops and executed branches are supported;
#' arbitrary compiled functions and matrix operations are not.
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
#' @section Higher-Order Derivatives:
#' \describe{
#'   \item{\code{\link{dual_variable_n}}}{Create a dual seeded for n-th order differentiation.}
#'   \item{\code{\link{deriv_n}}}{Extract the k-th derivative from a nested dual result.}
#'   \item{\code{\link{differentiate_n}}}{Compute f(x) and all derivatives up to order n.}
#' }
#'
#' @section Multi-Parameter Derivatives:
#' \describe{
#'   \item{\code{\link{D}}}{Composable total derivative operator. \code{D(f)}
#'     returns the derivative function; apply k times for k-th order tensors.}
#'   \item{\code{\link{gradient}}}{Gradient of a scalar-valued function.}
#'   \item{\code{\link{hessian}}}{Hessian matrix of a scalar-valued function.}
#'   \item{\code{\link{jacobian}}}{Jacobian matrix of a vector-valued function.}
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
#' @name nabla-package
#' @aliases nabla
"_PACKAGE"
