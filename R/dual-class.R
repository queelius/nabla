#' @title Dual Number Class for Automatic Differentiation
#' @description S4 class representing a dual number \eqn{a + b\varepsilon}
#'   where \eqn{\varepsilon^2 = 0}. The \code{value} slot holds the primal
#'   value and the \code{deriv} slot holds the tangent (derivative) component.
#'   Both slots accept \code{ANY} type to support nested duals for higher-order
#'   derivatives.
#' @slot value The primal (function) value. Numeric for first-order duals,
#'   or another dual for higher-order.
#' @slot deriv The tangent (derivative) component. Numeric for first-order
#'   duals, or another dual for higher-order.
#' @exportClass dual
setClass("dual", slots = list(value = "ANY", deriv = "ANY"))

# -- Constructors --------------------------------------------------------------

#' Create a dual number
#'
#' @param value The primal value (numeric or dual for nesting).
#' @param deriv The derivative component (numeric or dual for nesting).
#'   Defaults to 0.
#' @return A \code{dual} object.
#' @export
#' @examples
#' x <- dual(3, 1)
#' value(x)  # 3
#' deriv(x)  # 1
dual <- function(value, deriv = 0) {
  new("dual", value = value, deriv = deriv)
}

#' Create a dual variable (derivative seed = 1)
#'
#' Convenience constructor for the independent variable when computing
#' derivatives. Sets \code{deriv = 1} so that the output's derivative slot
#' contains \eqn{df/dx}.
#'
#' @param x A numeric value.
#' @return A \code{dual} with \code{value = x} and \code{deriv = 1}.
#' @export
#' @examples
#' x <- dual_variable(2)
#' deriv(x^2)  # 4 = 2*x evaluated at x=2
dual_variable <- function(x) {
  dual(x, 1)
}

#' Create a dual constant (derivative seed = 0)
#'
#' Wraps a numeric value as a dual with zero derivative, representing a
#' constant with respect to the differentiation variable.
#'
#' @param x A numeric value.
#' @return A \code{dual} with \code{value = x} and \code{deriv = 0}.
#' @export
#' @examples
#' k <- dual_constant(5)
#' deriv(k)  # 0
dual_constant <- function(x) {
  dual(x, 0)
}

# -- Generics and accessors ----------------------------------------------------

#' Extract the value (primal) part of a dual number
#'
#' @param d A \code{dual} object.
#' @return The value slot.
#' @export
#' @examples
#' value(dual(3, 1))  # 3
setGeneric("value", function(d) standardGeneric("value"))

#' @rdname value
#' @export
setMethod("value", "dual", function(d) d@value)

#' @rdname value
#' @export
setMethod("value", "numeric", function(d) d)

#' Extract the derivative (tangent) part of a dual number
#'
#' @param d A \code{dual} object.
#' @return The deriv slot.
#' @export
#' @examples
#' deriv(dual(3, 1))  # 1
setGeneric("deriv", function(d) standardGeneric("deriv"))

#' @rdname deriv
#' @export
setMethod("deriv", "dual", function(d) d@deriv)

#' @rdname deriv
#' @export
setMethod("deriv", "numeric", function(d) 0)

# -- Display -------------------------------------------------------------------

#' Display a dual number
#'
#' @param object A \code{dual} object.
#' @return Invisible \code{NULL}; called for side effect of printing.
#' @examples
#' x <- dual(3, 1)
#' x  # prints: <dual: 3 + 1*e>
#'
#' dv <- dual_vector(dual(1, 0), dual(2, 1))
#' dv  # prints: <dual_vector: 2 elements>
#'
#' @name dual-show
#' @aliases show,dual-method show,dual_vector-method
NULL

#' @rdname dual-show
#' @export
setMethod("show", "dual", function(object) {
  v <- format(value(object))
  d <- format(deriv(object))
  cat(sprintf("<dual: %s + %s*e>\n", v, d))
})

# -- Coercion -----------------------------------------------------------------

#' Coerce dual to numeric
#'
#' Extracts the primal value, discarding the derivative.
#'
#' @param x A \code{dual} object.
#' @param ... Ignored.
#' @return Numeric value.
#' @examples
#' x <- dual(3.14, 1)
#' as.numeric(x)  # 3.14
#'
#' @name dual-coerce
#' @aliases as.numeric,dual-method
NULL

#' @rdname dual-coerce
#' @export
setMethod("as.numeric", "dual", function(x, ...) {
  as.numeric(value(x))
})

# -- Predicate -----------------------------------------------------------------

#' Test whether an object is a dual number
#'
#' @param x An object.
#' @return Logical.
#' @examples
#' is_dual(dual(1, 0))  # TRUE
#' is_dual(42)           # FALSE
#' @export
is_dual <- function(x) {
  is(x, "dual")
}

# -- is.numeric: return TRUE so defensive checks pass -------------------------

#' Check if a dual number is numeric
#'
#' Returns \code{TRUE} for dual numbers so that defensive type checks pass.
#'
#' @param x A \code{dual} object.
#' @return \code{TRUE}.
#' @examples
#' is.numeric(dual(1, 0))  # TRUE
#'
#' @name dual-is-numeric
#' @aliases is.numeric,dual-method
NULL

#' @rdname dual-is-numeric
#' @export
setMethod("is.numeric", "dual", function(x) TRUE)

# -- c() method: collect duals into a dual_vector -----------------------------

#' Combine dual numbers into a dual_vector
#'
#' @param x A \code{dual} number.
#' @param ... Additional duals or numerics.
#' @param recursive Ignored.
#' @return A \code{dual_vector}.
#' @examples
#' x <- dual_variable(1)
#' y <- dual_variable(2)
#' dv <- c(x, y)
#' length(dv)  # 2
#'
#' @name dual-combine
#' @aliases c,dual-method
NULL

#' @rdname dual-combine
#' @export
setMethod("c", "dual", function(x, ..., recursive = FALSE) {
  args <- list(x, ...)
  args <- lapply(args, .as_dual)
  dual_vector(args)
})

# -- Internal helper: promote numeric to dual constant -------------------------

.as_dual <- function(x) {
  if (is(x, "dual")) x else dual(x, 0)
}

# -- Dual vector: a list-like container for multiple dual numbers ---------------

#' @title Dual Number Vector
#' @description A container for multiple dual numbers that supports indexing
#'   with \code{[} and \code{[[}, allowing log-likelihood functions to be
#'   written with \code{theta[1]}, \code{theta[2]} notation.
#' @slot .Data List of dual objects.
#' @exportClass dual_vector
setClass("dual_vector", contains = "list")

#' Create a vector of dual numbers
#'
#' Wraps a list of dual objects in a container that supports \code{[]} indexing
#' and \code{length()}, so that user functions can use natural
#' \code{theta[1]} notation.
#'
#' @param ... Dual objects, or a single list of dual objects.
#' @return A \code{dual_vector}.
#' @export
#' @examples
#' dv <- dual_vector(dual(1, 0), dual(2, 1))
#' length(dv)  # 2
#' value(dv[1])  # 1
dual_vector <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1L]]) && !is(args[[1L]], "dual")) {
    args <- args[[1L]]
  }
  new("dual_vector", args)
}

#' Indexing and length for dual_vector
#'
#' @param x A \code{dual_vector}.
#' @param i Numeric index.
#' @param j,drop,...  Ignored (present for generic compatibility).
#' @return A single \code{dual} for scalar index; a \code{dual_vector} for
#'   vector index; an integer for \code{length}.
#' @examples
#' dv <- dual_vector(dual(10, 1), dual(20, 0), dual(30, 0))
#' value(dv[1])   # 10
#' length(dv)     # 3
#'
#' @name dual_vector-access
#' @aliases [,dual_vector,numeric-method length,dual_vector-method
NULL

#' @rdname dual_vector-access
#' @export
setMethod("[", signature(x = "dual_vector", i = "numeric"),
  function(x, i, j, ..., drop = TRUE) {
    if (length(i) == 1L) {
      x[[i]]
    } else {
      dual_vector(methods::callNextMethod())
    }
  }
)

#' @rdname dual_vector-access
#' @export
setMethod("length", "dual_vector", function(x) {
  length(x@.Data)
})

#' @rdname dual-show
#' @export
setMethod("show", "dual_vector", function(object) {
  cat(sprintf("<dual_vector: %d elements>\n", length(object)))
  for (i in seq_along(object@.Data)) {
    cat(sprintf("  [%d] ", i))
    show(object[[i]])
  }
})
