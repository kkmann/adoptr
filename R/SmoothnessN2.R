#' Quadratic smoothness penalty term
#'
#' \code{SmoothnessN2} is a generic class for implementing a smoothness penalty
#' via the average squared first derivative of the stage two sample size function
#' \code{n2} of a two-stage design.
#' The only parameter is the width used for the finite differences, \code{h}.
#' The generic implementation only evluates \code{n2} in the interior of the
#' continuation region of a design.
#'
#' @slot distribution Data distribution
#' @slot prior Prior distribution
#' @slot h positive number giving the width of the central finite difference
#'     interval for approximating the first derivative.
#'
#' @exportClass SmoothnessN2
setClass("SmoothnessN2", representation(
        distribution = "DataDistribution",
        prior = "Prior",
        h = "numeric"
    ),
    contains = "ConditionalScore")



#' @param distribution see slot
#' @param prior see slot
#' @param h positive number, see slot \code{h}
#'
#' @rdname SmoothnessN2-class
#' @export
SmoothnessN2 <- function(distribution,
                         prior = ContinuousPrior(
                             pdf = function(x) stats::dunif(x, -5, 5),
                             support = c(-5, 5)),
                         h = .1)
    new("SmoothnessN2", distribution = distribution, prior = prior, h = h)




#' A generic implementation for arbitrary two-stage designs based on adaptive
#' Gaussian quadrature integration of the finite-differences approximation to
#' the first derivative is provided.
#' Custom subclasses of \code{Design} might implement this slightly different
#' (cf. [TODO: link to GQDesign implementation]).
#'
#' @param s an object of class \code{SmoothnessN2}
#' @param design the design to compute the smoothness term for
#' @param x1 first-stage outcome
#' @template dotdotdotTemplate
#'
#' @rdname SmoothnessN2-class
#' @export
setMethod("evaluate", signature("SmoothnessN2", "TwoStageDesign"),
          function(s, design, x1, ...) {
              res <- abs(n2(design, x1 + s@h) - n2(design, x1 - s@h))
              return(res / s@h)
          }
)

