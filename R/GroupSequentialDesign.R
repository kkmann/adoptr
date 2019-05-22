#' Group-sequential two-stage designs
#'
#' Group-sequential designs are a sub-class of the \code{TwoStageDesign}
#' class with constant stage-two sample size.
#' See \code{\link{TwoStageDesign}} for slot details.
#' Any group-sequential design can be converted to a fully flexible
#' \code{TwoStageDesign} (see examples section).
#'
#' @seealso \code{\link{TwoStageDesign}} for superclass and inherited methods
#'
#' @exportClass GroupSequentialDesign
setClass("GroupSequentialDesign",  contains = "TwoStageDesign")

#' @template c1f
#' @template c1e
#' @param n2_pivots numeric of length one, stage-two sample size
#' @param c2_pivots numeric vector, stage-two critical values on the integration
#' pivot points
#' @param order of the Gaussian uadrature rule to use for integration, set to
#' length(c2_pivots) if NULL, otherwise first value of c2_pivots is repeated
#' 'order'-times.
#' @template dotdotdot
#'
#' @examples
#' design <- GroupSequentialDesign(25, 0, 2, 25, c(1, 1.5, 2.5))
#' summary(design)
#'
#' @rdname GroupSequentialDesign-class
#' @export
GroupSequentialDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots, order = NULL, ...) {
    if (is.null(order)) {
        order <- length(c2_pivots)
        } else if (length(c2_pivots) != order) {
            c2_pivots <- rep(c2_pivots[1], order)
        }

     rule <- GaussLegendreRule(as.integer(order))

     tunable <- logical(8) # initialize to all false
     tunable[1:5] <- TRUE
     names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")

      new("GroupSequentialDesign", n1 = n1, c1f = c1f, c1e = c1e,
          n2_pivots = n2_pivots, c2_pivots = c2_pivots,
          x1_norm_pivots = rule$nodes, weights = rule$weights, tunable = tunable)
}





#' @rdname n
#' @export
setMethod("n2", signature("GroupSequentialDesign", "numeric"),
          function(d, x1, round = TRUE, ...) {
              n2 <- ifelse(x1 < d@c1f | x1 > d@c1e, 0, d@n2_pivots)
              if (round)
                  n2 <- round(n2)
              return(n2)
          }
)

#' @param n1 stage one sample size or \code{GroupSequentialDesign} object to convert
#'   (overloaded from \code{\link{TwoStageDesign}})
#'
#' @examples
#' TwoStageDesign(design)
#'
#' @rdname GroupSequentialDesign-class
#' @export
setMethod("TwoStageDesign", signature("GroupSequentialDesign"),
     function(n1, ...){
         tunable <- logical(8) # initialize to all false
         tunable[1:5] <- TRUE
         names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")
         new("TwoStageDesign", n1 = n1@n1, c1f = n1@c1f, c1e = n1@c1e,
                        n2_pivots = rep(n1@n2_pivots, length(n1@weights)),
                        c2_pivots = n1@c2_pivots,
                        x1_norm_pivots = n1@x1_norm_pivots, weights = n1@weights,
                        tunable = tunable)
})
