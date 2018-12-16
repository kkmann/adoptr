#' Design with Gaussian quadrature integration rule
#'
#'     gq_design simply creates a TwoStageDesign object with a Gaussian
#'     quadrature numerical integration rule. If length(n2) = 1 and
#'     length(c2) > 1 an object of class \code{GSDesign} is created.
#'
#' @param n1 stage-one sample size
#' @param c1f early stopping for futility boundary
#' @param c1e early stopping for efficacy boundary
#' @param n2 second-stage sample size. vector for adaptive and
#'           integer for group-sequential design
#' @param c2 secnod-stage decision boundary
#' @param order order (i.e. number of pivot points in the interior of [c1f, c1e])
#'     of the Gaussian quadrature rule to use for integration
#' @param ... futher arguments
#'
#' @export
setGeneric("gq_design", function(n1, c1f, c1e, n2, c2, order, ...)
    standardGeneric("gq_design"))

#' @describeIn gq_design Method to create design with Gaussian quadrature rule
#'             and Legendre nodes
#' @export
setMethod("gq_design", signature("numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric"),
          function(n1, c1f, c1e, n2, c2, order, ...) {
              rule <- GaussLegendreRule(order)
              if(length(n2) == length(c2) & length(n2) == order) {
                  TwoStageDesign(n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2,
                                 c2_pivots = c2, x1_norm_pivots = rule$nodes,
                                 weights = rule$weights)
              } else if (length(n2) == 1 & length(c2) > 1 & length(c2) == order) {
                  GSDesign(n1 = n1, c1f = c1f, c1e = c1e,
                           n2_continue = n2, c2_pivots = c2,
                           x1_norm_pivots = rule$nodes, weights = rule$weights)
              } else {
                  stop("parameter lengths do not fit!")
              }
          }
)

