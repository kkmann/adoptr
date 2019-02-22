#' Group-sequential design
#'
#' The class \code{GroupSequentialDesign} allows the implementation
#' and usage of group-sequential designs. These are two-stage designs
#' with a constant stage-two sample size function.
#' Therefore, the length of \code{n2_pivots} is required to be 1 in this
#' class.
#' \code{GroupSequentialDesign} is a subclass of \link{TwoStageDesign}.
#'
#' @slot n1 stage-one sample size
#' @slot c1f early stopping for futility boundary
#' @slot c1e early stopping for efficacy boundary
#' @slot n2_pivots stage-two sample size upon continuation
#' @slot c2_pivots vector of length order giving the values of c2 at the
#'     pivot points of the numeric integration rule [TODO: these are not available during construction]
#' @slot x1_norm_pivots normalized pivots for integration rule (in [-1, 1])
#' @slot weights weights of conditional score values at x1_norm_pivots for
#'     approximating the integral over x1.
#' @slot tunable named logical vector indicating whether corresponding slot is considered a tunable parameter
#'
#' @exportClass GroupSequentialDesign
setClass("GroupSequentialDesign",  contains = "TwoStageDesign")





#' Create an object of class \code{GroupSequentialDesign}
#'
#' If \code{order} is not defined, the length of the parameter \code{c2_pivots}
#' defines the order of the Gaussian quadrature integration rule.
#'
#' @param n1 cf. slot
#' @param c1f cf. slot
#' @param c1e cf. slot
#' @param n2_pivots stage-two sample size
#' @param c2_pivots stage-two rejection boundary
#' @param order order of Gaussian quadrature
#' @param ... further optional arguments
#'
#' @return an object of class \code{GroupSequentialDesign}
#'
#' @rdname GroupSequentialDesign-class
#' @export
GroupSequentialDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots, order = NULL, ...) {
    if(is.null(order)) {
        order <- length(c2_pivots)
        } else{
            c2_pivots <- rep(c2_pivots[1], order)
        }

     rule <- GaussLegendreRule(order)

     tunable <- logical(8) # initialize to all false
     tunable[1:5] <- TRUE
     names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")

      new("GroupSequentialDesign", n1 = n1, c1f = c1f, c1e = c1e,
          n2_pivots = n2_pivots, c2_pivots = c2_pivots,
          x1_norm_pivots = rule$nodes, weights = rule$weights, tunable = tunable)
}





#' @param x1 stage-one outcome
#' @param d object of class \code{GSDesign}
#' @param round logical, should integer sample size or real sample size be
#'    returned?
#'
#' @rdname GroupSequentialDesign-class
#' @export
setMethod("n2", signature("GroupSequentialDesign", "numeric"),
          function(d, x1, round = TRUE, ...) {
              n2 <- ifelse(x1 < d@c1f | x1 > d@c1e, 0, d@n2_pivots)
              if (round)
                  n2 <- round(n2)
              return(n2)
          }
)



#' Convert a group-sequential design to a two-stage design
#'
#' @param tunable c.f. slot
#'
#' @rdname GroupSequentialDesign-class
#' @export
setMethod("TwoStageDesign", signature("GroupSequentialDesign"),
     function(d, tunable, ...){
         tunable <- logical(8) # initialize to all false
         tunable[1:5] <- TRUE
         names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")
         new("TwoStageDesign", n1 = d@n1, c1f = d@c1f, c1e = d@c1e,
                        n2_pivots = rep(d@n2_pivots, length(d@weights)),
                        c2_pivots = d@c2_pivots,
                        x1_norm_pivots = d@x1_norm_pivots, weights = d@weights,
                        tunable = tunable)
})

