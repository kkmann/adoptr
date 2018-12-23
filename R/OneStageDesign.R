#' One-stage design
#'
#' [ToDo]
#'
#' @slot n1 sample size
#' @slot c1f rejection boundary
#' @slot c1e rejection boundary
#' @slot n2_pivots constant 0
#' @slot c2_pivots \code{+/- Inf}
#' @slot x1_norm_pivots is ignored for OneStageDesign
#' @slot weights is ignored for OneStageDesign
#'
#' @exportClass OneStageDesign
setClass("OneStageDesign",  contains = "TwoStageDesign")



#' @param n sample size (stage-one sample size)
#' @param c rejection boundary (c = c1f = c1e)
#'
#' @rdname OneStageDesign-class
#' @export
OneStageDesign <- function(n, c) {
    new("OneStageDesign", n1 = n, c1f = c, c1e = c, n2_pivots = 0,
    c2_pivots = NaN, x1_norm_pivots = NaN, weights = NaN)
}




#' @param x object to get parameters from
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("tunable_parameters", signature("OneStageDesign"),
          function(x, ...) c(x@n1, x@c1f))




#' @param params vector of design parameters (must be in same order as returned
#'     by \code{as.numeric(OneStageDesign)})
#' @param object object to update
#' @param ... further optional arguments
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("update", signature("OneStageDesign"),
    function(object, params, ...) {
        if(length(params) != 2)
            stop("parameter length does not fit")
        new("OneStageDesign",
            n1  = params[1],
            c1f = params[2],
            c1e = params[2],
            n2_pivots    = 0,
            c2_pivots      = NaN,
            x1_norm_pivots = NaN,
            weights = NaN)
    })



#' @param x1 stage-one outcome
#' @param d object of class \code{OneStageDesign}
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("n2", signature("OneStageDesign", "numeric"),
          function(d, x1, ...) 0 )


#' @rdname OneStageDesign-class
#' @export
setMethod("c2", signature("OneStageDesign", "numeric"),
          function(d, x1, ...) ifelse(x1 <= d@c1f, Inf, -Inf) )





#' Convert a one-stage design to a two-stage design
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("TwoStageDesign", signature("OneStageDesign"),
     function(d, ...){
        new("TwoStageDesign", n1 = d@n1, c1f = d@c1f, c1e = d@c1f,
                        n2_pivots = rep(0, length(d@weights)),
                        c2_pivots = rep(Inf, length(d@weights)),
                        x1_norm_pivots = d@x1_norm_pivots, weights = d@weights)
})

