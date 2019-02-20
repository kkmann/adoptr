#' One-stage design
#'
#' \code{OneStageDesign} allows the usage of a single-stage design.
#' In this case without a second stage, there only exist the first-stage sample
#' size \code{n1} and the first-stage rejection boundary \code{c1f}.
#' No other parameters have to be set.
#' \code{OneStageDesign} is a subclass of \link{TwoStageDesign}.
#'
#' @slot n1 sample size
#' @slot c1f rejection boundary
#' @slot c1e rejection boundary
#' @slot n2_pivots constant 0
#' @slot c2_pivots \code{+/- Inf}
#' @slot x1_norm_pivots is ignored for OneStageDesign
#' @slot weights is ignored for OneStageDesign
#' @slot tunable is ignored for OneStageDesign
#' @slot rounded logical that indicates whether rounded n1-value should be used
#'
#' @exportClass OneStageDesign
setClass("OneStageDesign",  contains = "TwoStageDesign")



#' @param n sample size (stage-one sample size)
#' @param c rejection boundary (c = c1f = c1e)
#'
#' @rdname OneStageDesign-class
#' @export
OneStageDesign <- function(n, c) {
    tunable <- logical(8)
    tunable[1:2] <- TRUE
    names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")
    new("OneStageDesign", n1 = n, c1f = c, c1e = c, n2_pivots = 0,
    c2_pivots = NaN, x1_norm_pivots = NaN, weights = NaN,
    tunable = tunable, rounded = FALSE)
}



setMethod("update", signature("TwoStageDesign"),
          function(object, params, ...) {
              tunable_names <- names(object@tunable)[object@tunable]
              res <- object
              idx <- 1
              for (i in 1:length(tunable_names)) {
                  slotname <- tunable_names[i]
                  k <- length(slot(object, name = slotname))
                  slot(res, name = slotname) <- params[idx:(idx + k - 1)]
                  idx <- idx + k
              }
              return(res)
          })



#' @param params vector of design parameters (must be in same order as returned
#'     by \code{as.numeric(OneStageDesign)})
#' @param object object to update
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("update", signature("OneStageDesign"),
          function(object, params, ...) {
              tunable_names <- names(object@tunable)[object@tunable]
              res <- object
              idx <- 1
              for (i in 1:length(tunable_names)) {
                  slotname <- tunable_names[i]
                  k <- length(slot(object, name = slotname))
                  slot(res, name = slotname) <- params[idx:(idx + k - 1)]
                  idx <- idx + k
              }
              res@c1e <- res@c1f
              return(res)
          })



#' @param x1 stage-one outcome
#' @param d object of class \code{OneStageDesign}
#' @param ... further optional arguments
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
#' @param rounded logical that indicates whether rounded n1-value should be used
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("TwoStageDesign", signature("OneStageDesign"),
     function(d, rounded = FALSE, ...){
         tunable <- rep(TRUE, 2)
         names(tunable) <- c("n1", "c1f")
         new("TwoStageDesign",
             n1 = d@n1,
             c1f = d@c1f - .01, # needs to be done for interpolation
             c1e = d@c1f + .01, # needs to be done for interpolation
             n2_pivots = rep(0, 2),
             c2_pivots = c(10, -10),
             x1_norm_pivots = c(-.5, .5),
             weights = c(1, 1),
             tunable = tunable,
             rounded = rounded)
})



#' plot() is not defined for one stage designs - what should be plotted?
#'
#' @param x not used
#' @param y not used
#' @param k number of points to use for plotting
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("plot", signature("OneStageDesign"),
          function(x, y = NULL, rounded = TRUE, ..., k = 100)
              stop("plot method is only defined for two-stage designs!")
          )
