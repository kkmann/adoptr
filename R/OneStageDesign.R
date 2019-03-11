#' One-stage designs
#'
#' \code{OneStageDesign} implements a one-stage design as special case of
#' a two-stage design, i.e. as sub-class of \code{\link{TwoStageDesign}}.
#' This is possible by defining n2 = 0, c:=c1f=c1e, c2(x1) = ifelse(x1 < c, Inf, -Inf).
#' No integration pivots etc are required (set to NaN).
#' Note that the default \code{\link[plot,TwoStageDesign-method]{plot}} method
#' is not supported for \code{OneStageDesign} objects.
#'
#' @seealso \code{\link{TwoStageDesign}}, \code{\link{GroupSequentialDesign}}
#'
#' @exportClass OneStageDesign
setClass("OneStageDesign",  contains = "TwoStageDesign")

#' @param n sample size (stage-one sample size)
#' @param c rejection boundary (c = c1f = c1e)
#'
#' @examples
#' design <- OneStageDesign(30, 1.96)
#' summary(design)
#' design <- TwoStageDesign(design)
#' summary(design)
#'
#' @rdname OneStageDesign-class
#' @export
OneStageDesign <- function(n, c) {
    tunable <- logical(8)
    tunable[1:2] <- TRUE
    names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")
    new("OneStageDesign", n1 = n, c1f = c, c1e = c, n2_pivots = 0,
    c2_pivots = NaN, x1_norm_pivots = NaN, weights = NaN,
    tunable = tunable)
}





#' @rdname tunable_parameters
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





#' @rdname n
#' @export
setMethod("n2", signature("OneStageDesign", "numeric"),
          function(d, x1, ...) 0 )





#' @rdname critical-values
#' @export
setMethod("c2", signature("OneStageDesign", "numeric"),
          function(d, x1, ...) ifelse(x1 <= d@c1f, Inf, -Inf) )





#' @rdname OneStageDesign-class
#' @export
setMethod("TwoStageDesign", signature("OneStageDesign"),
     function(d, ...){
         tunable <- rep(TRUE, 2)
         names(tunable) <- c("n1", "c1f")
         new("TwoStageDesign",
             n1  = d@n1,
             c1f = d@c1f - .01, # needs to be done for interpolation
             c1e = d@c1f + .01, # needs to be done for interpolation
             n2_pivots = rep(0, 2),
             c2_pivots = c(10, -10),
             x1_norm_pivots = c(-.5, .5),
             weights = c(1, 1),
             tunable = tunable)
})



#' plot() is not defined for one stage designs
#'
#' @param x not used
#' @param y not used
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("plot", signature("OneStageDesign"),
          function(x, y, ...)
              stop("plot method is only defined for two-stage designs!")
          )
