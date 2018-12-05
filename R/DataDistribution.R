setClass("DataDistribution")

setGeneric("probability_density_function", function(dist, x, n, theta, ...) standardGeneric("probability_density_function"))

setGeneric("cumulative_distribution_function", function(dist, x, n, theta, ...) standardGeneric("cumulative_distribution_function"))

setMethod("quantile", signature("DataDistribution"), function(x, probs, n, theta, ...) stop("not implemented"))








setClass("Normal", representation(
    dummy = "logical" # needed to make this non-abstract
    ),
    contains = "DataDistribution")

Normal <- function() new("Normal", dummy = FALSE)

setMethod("probability_density_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) stats::dnorm(x, mean = sqrt(n) * theta, sd = 1) )

setMethod("cumulative_distribution_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) stats::pnorm(x, mean = sqrt(n) * theta, sd = 1) )

setMethod("quantile", signature("DataDistribution"),
          function(x, probs, n, theta, ...) stats::qnorm(probs, mean = sqrt(n) * theta, sd = 1) )
