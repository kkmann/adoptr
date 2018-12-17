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
    new("GSDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
        c2_pivots = c2_pivots, x1_norm_pivots = x1_norm_pivots, weights = weights)
}





#' @param x object to get parameters from
#'
#' @rdname GSDesign-class
#' @export
setMethod("tunable_parameters", signature("GSDesign"),
          function(x, ...) c(x@n1, x@c1f, x@c1e, x@n2_pivots, x@c2_pivots))




#' @param params vector of design parameters (must be in same order as returned
#'     by \code{as.numeric(design)})
#' @param object object to update
#'
#' @rdname GSDesign-class
#' @export
setMethod("update", signature("GSDesign"),
    function(object, params, ...) {
        k <- length(object@weights)
        if( (length(params) - 4) != k)
            stop("parameter length does not fit")
        new("GSDesign",
            n1  = params[1],
            c1f = params[2],
            c1e = params[3],
            n2_pivots    = params[4],
            c2_pivots      = params[5:(length(params))],
            x1_norm_pivots = object@x1_norm_pivots,
            weights = object@weights)
    })



#' @param x1 stage-one outcome
#' @param d object of class \code{GSDesign}
#'
#' @rdname GSDesign-class
#' @export
setMethod("n2", signature("GSDesign", "numeric"),
          function(d, x1, ...) ifelse(x1 < d@c1f | x1 > d@c1e, 0, d@n2_pivots) )


#' Return smootheness of a group-sequential design as 0.
#'
#' @param s an object of class \code{SmoothnessN2}
#' @param design an object of class \code{GSDesign}
#'
#' @rdname GSDesign-class
#' @export
setMethod("evaluate", signature("SmoothnessN2", "GSDesign"),
          function(s, design, ...) 0 )



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

