#' @exportClass ConditionalSampleSize
setClass("ConditionalSampleSize", contains = "ConditionalScore")

#' @export
ConditionalSampleSize <- function(dist, prior) new("ConditionalSampleSize", distribution = dist, prior = prior)

setMethod("evaluate", signature("ConditionalSampleSize", "Design"),
          function(s, design, x1, ...) n(design, x1) )
