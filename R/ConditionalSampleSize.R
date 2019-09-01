#' (Conditional) Sample Size of a Design
#'
#' This score simply evaluates \code{n(d, x1)} for a design \code{d} and the
#' first-stage outcome \code{x1}.
#' The data distribution and prior are only relevant when it is integrated.
#'
#' @template dist
#' @template prior
#' @template design
#' @template s
#' @template x1
#' @template optimization
#' @template dotdotdot
#'
#' @seealso \link{Scores}
#'
#' @examples
#' design <- TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L)
#' prior  <- PointMassPrior(.4, 1)
#'
#' css   <- ConditionalSampleSize()
#' evaluate(css, design, c(0, .5, 3))
#'
#' ess   <- ExpectedSampleSize(Normal(), prior)
#'
#' # those two are equivalent
#' evaluate(ess, design)
#' evaluate(expected(css, Normal(), prior), design)
#'
#' @aliases ConditionalSampleSize
#' @exportClass ConditionalSampleSize
setClass("ConditionalSampleSize", contains = "ConditionalScore")



#' @template label
#'
#' @rdname ConditionalSampleSize-class
#' @export
ConditionalSampleSize <- function(label = "n(x1)") new("ConditionalSampleSize", label = label)



#' @rdname ConditionalSampleSize-class
#' @export
ExpectedSampleSize <- function(dist, prior, label = "E[n(x1)]") {
    expected(ConditionalSampleSize(), dist, prior, label = label)
}



#' @rdname ConditionalSampleSize-class
#' @export
setMethod("evaluate", signature("ConditionalSampleSize", "TwoStageDesign"),
          function(s, design, x1, optimization = FALSE, ...) {
              n(design, x1, round = !optimization)
          })
