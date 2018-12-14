#' Two-stage designs
#'
#' [ToDo]
#'
#' @slot n1 stage-one sample size
#' @slot c1f early stopping for futility boundary
#' @slot c1e early stopping for efficacy boundary
#' @slot n2_pivots vector of length order giving the values of n2 at the
#'     pivot points of the numeric integration rule [TODO: these are not available during construction]
#' @slot c2_pivots vector of length order giving the values of c2 at the
#'     pivot points of the numeric integration rule [TODO: these are not available during construction]
#' @slot x1_norm_pivots normalized pivots for integration rule (in [-1, 1])
#' @slot weights weights of conditional score values at x1_norm_pivots for
#'     approximating the integral over x1.
#'
#' @exportClass TwoStageDesign
setClass("TwoStageDesign", representation(
        n1        = "numeric",
        c1f       = "numeric",
        c1e       = "numeric",
        n2_pivots = "numeric",
        c2_pivots = "numeric",
        x1_norm_pivots = "numeric",
        weights   = "numeric"
    ))



#' @param ... optional arguments depending on implementation
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("TwoStageDesign", function(...) standardGeneric("TwoStageDesign"))

#' @param n1 cf. slot
#' @param c1f cf. slot
#' @param c1e cf. slot
#' @param n2_pivots cf. slot
#' @param c2_pivots cf. slot
#' @param x1_norm_pivots cf. slot
#' @param weights cf. slot
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("TwoStageDesign", signature("numeric"),
     function(n1, c1f, c1e, n2_pivots, c2_pivots, x1_norm_pivots, weights) {
        if (any(diff(sapply(list(n2_pivots, c2_pivots, x1_norm_pivots, weights), length)) != 0))
            stop("pivots and weights must all be of the same length")
        if (any(x1_norm_pivots < -1) | any(x1_norm_pivots > 1))
            stop("x1_norm_pivots must be in [-1, 1], is scaled automatically")
        if (any(weights <= 0))
            stop("weights must be positive")
        new("TwoStageDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
            c2_pivots = c2_pivots, x1_norm_pivots = x1_norm_pivots, weights = weights)
})



#' @param x object to get parameters from
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("tunable_parameters", function(x, ...) standardGeneric("tunable_parameters"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("tunable_parameters", signature("TwoStageDesign"),
          function(x, ...) c(x@n1, x@c1f, x@c1e, x@n2_pivots, x@c2_pivots))





#' @param params vector of design parameters (must be in same order as returned
#'     by \code{as.numeric(design)})
#' @param object object to update
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("update", signature("TwoStageDesign"),
    function(object, params, ...) {
        k <- length(object@weights)
        if( ((length(params) - 3) / 2) != k)
            stop("parameter length does not fit")
        new("TwoStageDesign",
            n1  = params[1],
            c1f = params[2],
            c1e = params[3],
            n2_pivots = params[4:(3 + k)],
            c2_pivots = params[(4 + k):(length(params))],
            x1_norm_pivots = object@x1_norm_pivots,
            weights = object@weights)
    })


#' @param x1 stage-one outcome
#' @param d design object
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("n2", function(d, x1, ...) standardGeneric("n2"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("n2", signature("TwoStageDesign", "numeric"),
          function(d, x1, ...) ifelse(x1 < d@c1f | x1 > d@c1e, 0, 1) *
              pmax(0, stats::approx(scaled_integration_pivots(d), d@n2_pivots, xout = x1, method = "linear", rule = 2)$y) )



#' @rdname TwoStageDesign-class
#' @export
setGeneric("n", function(d, x1, ...) standardGeneric("n"))

#' @describeIn TwoStageDesign overall sample size given stage-one outcome
#' @export
setMethod("n", signature("TwoStageDesign", "numeric"), function(d, x1, ...) n2(d, x1, ...) + d@n1 )



#' @rdname TwoStageDesign-class
#' @export
setGeneric("c2", function(d, x1, ...) standardGeneric("c2"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("c2", signature("TwoStageDesign", "numeric"),
    function(d, x1, ...) stats::approx(scaled_integration_pivots(d), d@c2_pivots, xout = x1, method = "linear", rule = 2)$y *
        ifelse(x1 < d@c1f, Inf, 1) * ifelse(x1 > d@c1e, -Inf, 1) )



#' @rdname TwoStageDesign-class
#' @export
setGeneric("scaled_integration_pivots", function(d, ...) standardGeneric("scaled_integration_pivots"))

#' @describeIn TwoStageDesign get the actual pivots points (knots) of the numerical integration routine
#'     rule.
setMethod("scaled_integration_pivots", signature("TwoStageDesign"),
          function(d, ...){
              h <- (d@c1e - d@c1f) / 2
              return(h * d@x1_norm_pivots + (h + d@c1f))
          })



#' Show method for TwoStageDesign objects
#'
#' Only states the class itself. More information is given by the respective
#' \code{summary} method, cf. \link{TwoStageDesign-class}.
#'
#' @param object design to show
#'
#' @export
setMethod("show", signature(object = "TwoStageDesign"),
          function(object) cat("TwoStageDesign"))




#' Plot TwoStageDesign with optional set of conditional scores
#'
#' [TODO]
#'
#' @param x design to plot
#' @param y not used
#' @param k number of points to use for plotting
#' @param ... optional named (Un)ConditionalScores to plot / include in the summary
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("plot", signature(x = "TwoStageDesign"),
          function(x, y = NULL, ..., k = 100) {
              scores <- list(...)
              if (!all(sapply(scores, function(s) is(s, "ConditionalScore"))))
                  stop("optional arguments must be ConditionalScores")
              opts  <- graphics::par(mfrow = c(1, length(scores) + 2))
              x1    <- seq(x@c1f - (x@c1e - x@c1f)/5, x@c1e + (x@c1e - x@c1f)/5, length.out = k)
              plot(x1, n(x, x1), 'l', ylim = c(0, 1.05 * max(n(x, x1))),
                   main = "Overall sample size", ylab = "")
              plot(x1, c2(x, x1), 'l', main = "Stage-two critical value", ylab = "")
              if (length(scores) > 0) {
                  for (i in 1:length(scores)) {
                      plot(x1, evaluate(scores[[i]], x, x1), 'l', main = names(scores[i]),
                           ylab = "")
                  }
              }
              graphics::par(opts)
          })



#' Summarize TwoStageDesign objects with optional set of scores
#'
#' [TODO]
#'
#' @param object design object to plot
#' @param ... optinal additional named UnconditionalScores
#'
#' @export
setMethod("summary", signature("TwoStageDesign"),
          function(object, ...) {
              scores <- list(...)
              if (!all(sapply(scores, function(s) is(s, "UnconditionalScore"))))
                  stop("optional arguments must be UnconditionalScores")
              res <- list(
                  design = object,
                  scores = sapply(scores, function(s) evaluate(s, object))
              )
              names(res$scores) <- names(scores)
              class(res) <- c("TwoStageDesignSummary", "list")
              return(res)
          })


#' Print obejct of class TwoStageDesignSummary
#'
#' @param x object to print
#' @param ... unused
#'
#' @export
print.TwoStageDesignSummary <- function(x, ...) {
    cat("TwoStageDesign with:\n\r")
    cat(sprintf("     n1: %6.1f\n\r", x$design@n1))
    cat(sprintf("    c1f: %6.1f\n\r", x$design@c1f))
    cat(sprintf("    c1e: %6.1f\n\r", x$design@c1e))
    if (length(x$scores) > 0) {
        cat("Unconditional scores:\n\r")
        for (i in 1:length(x$scores)) {
            cat(sprintf("    %s: %7.2f\n\r", names(x$scores)[i], x$scores[i]))
        }
        cat("\n\r")
    }
}
