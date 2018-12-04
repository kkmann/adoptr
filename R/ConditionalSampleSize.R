setClass("ConditionalSampleSize", contains = "ConditionalScore")

ConditionalSampleSize <- function(prior) new("ConditionalSampleSize", prior = prior)

setMethod("evaluate", signature("ConditionalSampleSize", "Design"),
          function(s, design, x1, ...) n(design, x1) )
