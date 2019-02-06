#' Group-sequential design
#'
#' [ToDo]
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
#'
#' @exportClass GSDesign
setClass("GSDesign",  contains = "TwoStageDesign")



#' @param n1 cf. slot
#' @param c1f cf. slot
#' @param c1e cf. slot
#' @param n2_pivots cf. slot
#' @param c2_pivots cf. slot
#' @param x1_norm_pivots cf. slot
#' @param weights cf. slot
#' @param ... further optional arguments
#'
#' @rdname GSDesign-class
#' @export
GSDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots, x1_norm_pivots, weights) {
    if (any(diff(sapply(list(c2_pivots, x1_norm_pivots, weights), length)) != 0))
        stop("pivots and weights must all be of the same length")
    if (any(x1_norm_pivots < -1) | any(x1_norm_pivots > 1))
        stop("x1_norm_pivots must be in [-1, 1], is scaled automatically")
    if (any(weights <= 0))
        stop("weights must be positive")
    tunable <- logical(8) # initialize to all false
    tunable[1:5] <- TRUE
    names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")
        new("GSDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
        c2_pivots = c2_pivots, x1_norm_pivots = x1_norm_pivots, weights = weights,
        tunable = tunable, rounded = FALSE)
}




#' @param x1 stage-one outcome
#' @param d object of class \code{GSDesign}
#'
#' @rdname GSDesign-class
#' @export
setMethod("n2", signature("GSDesign", "numeric"),
          function(d, x1, ...) {
              res <- ifelse(x1 < d@c1f | x1 > d@c1e, 0, d@n2_pivots)
              if(d@rounded)
                  res <- round(res)
              return(res)
          }
)



#' Convert a group-sequential design to a two-stage design
#'
#' @rdname GSDesign-class
#' @export
setMethod("TwoStageDesign", signature("GSDesign"),
     function(d, ...){
        new("TwoStageDesign", n1 = d@n1, c1f = d@c1f, c1e = d@c1e,
                        n2_pivots = rep(d@n2_pivots, length(d@weights)),
                        c2_pivots = d@c2_pivots,
                        x1_norm_pivots = d@x1_norm_pivots, weights = d@weights)
})

