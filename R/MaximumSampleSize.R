#' Maximum Sample Size of a Design
#'
#' This score evaluates \code{max(n(d))} for a design \code{d}.
#'
#' @template s
#' @template design
#' @template optimization
#' @template dotdotdot
#'
#' @seealso \link{Scores} for general scores and \link{ConditionalSampleSize}
#'  for evaluating the sample size point-wise.
#'
#' @examples
#' design <- TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L)
#' mss    <- MaximumSampleSize()
#' evaluate(mss, design)
#'
#' @aliases MaximumSampleSize
#' @exportClass MaximumSampleSize
setClass("MaximumSampleSize", contains = "UnconditionalScore")



#' @template label
#'
#' @rdname MaximumSampleSize-class
#' @export
MaximumSampleSize <- function(label = "max(n(x1))") new("MaximumSampleSize", label = label)



#' @rdname MaximumSampleSize-class
#' @export
setMethod("evaluate",
          signature("MaximumSampleSize", "TwoStageDesign"),
          function(s, design, optimization = FALSE, ...) {
            x1 <- seq(design@c1f, design@c1e, length.out = 1000)
            ss <- sapply(x1, function(z) adoptr:::n(design, z, round = !optimization))
            return(max(ss))
          })
