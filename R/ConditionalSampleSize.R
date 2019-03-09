#' Conditional sample size of a design given stage-one outcome
#'
#' This score simply evaluates \code{n(d, x1)}.
#' The data distribution and prior are only relevant when it is integrated.
#'
#' @template dist
#' @template prior
#'
#' @seealso The method \code{\link{evaluate}} provides evaluation of the
#'    \code{ConditionalSampleSize}.
#'
#' @aliases ConditionalSampleSize
#' @exportClass ConditionalSampleSize
setClass("ConditionalSampleSize", contains = "ConditionalScore")

#' @examples
#' css <- ConditionalSampleSize(dist = Normal(), prior = PointMassPrior(.4, 1))
#'
#' @rdname ConditionalSampleSize-class
#' @export
ConditionalSampleSize <- function(dist, prior) new("ConditionalSampleSize", distribution = dist, prior = prior)


#' @examples
#' design <- TwoStageDesign(50, .0, 2.0, 50, 2.0) # dummy design
#'  evaluate(
#'     ConditionalSampleSize(Normal(), PointMassPrior(.3, 1)),
#'     design,
#'     x1 = 3
#'  ) # 50
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("ConditionalSampleSize", "TwoStageDesign"),
          function(s, design, x1, optimization = FALSE, ...) {
              n(design, x1, round = !optimization)
          })
