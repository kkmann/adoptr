#' Quadratic smoothness penalty term
#'
#' \code{Smoothness_n2} is a generic class for implementing a smoothness penalty
#' via the average squared first derivative of the stage two sample size function
#' \code{n2} of a two-stage design.
#' The only parameter is the width used for the finite differences, \code{h}.
#' The generic implementation only evluates \code{n2} in the interior of the
#' continuation region of a design.
#'
#' @slot h positive number giving the width of the central finite difference
#'     interval for approximating the first derivative.
#'
#' @exportClass Smoothness_n2
setClass("Smoothness_n2", representation(
        h = "numeric"
    ),
    contains = "IntegralScore")



#' @param h positive number, see slot \code{h}
#'
#' @rdname Smoothness_n2-class
#' @export
Smoothness_n2 <- function(h = sqrt(.Machine$double.eps)) new("Smoothness_n2", h = h)



#' A generic implementation for arbitrary two-stage designs based on adaptive
#' Gaussian quadrature integration of the finite-differences approximation to
#' the first derivative is provided.
#' Custom subclasses of \code{Design} might implement this slightly different
#' (cf. [TODO: link to GQDesign implementation]).
#'
#' @param s an object of class \code{Smoothness_n2}
#' @param design the design to compute the smoothness term for
#' @param specific logical, should a design-specific implementation be used?
#'     defaults to \code{TRUE}. If \code{TRUE}, looks for design-specific
#'     \code{.evaluate} method and calls it
#' @template dotdotdotTemplate
#'
#' @rdname Smoothness_n2-class
#' @export
setMethod("evaluate", signature("Smoothness_n2", "Design"),
          function(s, design, specific = TRUE, ...) {
              if (specific) {
                  # use design-specific implementation
                  return(.evaluate(s, design, ...))
              } else {
                  # use generic approach
                  # integrand is the finite difference approximation of the
                  # squared derivative
                  integrand <- function(x1) ((n2(design, x1 + s@h/2) - n2(design, x1 - s@h/2))/s@h)^2
                  x1_bounds <- early_stopping_bounds(design, ...) + c(s@h, -s@h)
                  # use adaptive quadrature to integrate - only relies on generic interface
                  # provided by 'Design', no special optimization for particular
                  # design implementation
                  return(1 / diff(x1_bounds) * stats::integrate(integrand, x1_bounds[1], x1_bounds[2])$value)
              }
          })
