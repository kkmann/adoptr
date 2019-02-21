#' Conditional sample size of a design given stage-one outcome
#'
#' This score simply evaluates \code{n(d, x1)}.
#' The data distribution and prior are only relevant when it is integrated.
#'
#' @details See documentation in \link{ConditionalScore-class}
#'     for further details and inherited methods.
#'
#' @template ConditionalScoreTemplate
#' @param dist data distribution
#' @param prior prior distribution
#'
#' @exportClass ConditionalSampleSize
setClass("ConditionalSampleSize", contains = "ConditionalScore")

#' @describeIn ConditionalSampleSize-class constructor
#' @export
ConditionalSampleSize <- function(dist, prior) new("ConditionalSampleSize", distribution = dist, prior = prior)

#' @param optimization logical, if TRUE uses a relaxation to real parameters of
#'    the underlying design; used for smooth optimization.
#' @rdname ConditionalSampleSize-class
setMethod("evaluate", signature("ConditionalSampleSize", "TwoStageDesign"),
          function(s, design, x1, optimization = FALSE, ...) {
              n(design, x1, round = !optimization)
          })
