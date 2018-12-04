setClass("DataDistribution")

setGeneric("probability_density_function", function(dist, x1, n1, theta, ...) standardGeneric("probability_density_function"))

setGeneric("cumulative_distribution_function", function(dist, x1, n1, theta, ...) standardGeneric("cumulative_distribution_function"))

setMethod("quantile", signature("DataDistribution"), function(x, probs, n1, theta, ...) stop("not implemented"))


setClass("Normal", representation(
    dummy = "logical" # needed to make this non-abstract
    ),
    contains = "DataDistribution")

Normal <- function() new("Normal", dummy = FALSE)

setMethod("probability_density_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x1, n1, theta, ...) dnorm(x1, mean = sqrt(n1) * theta, sd = 1) )

setMethod("cumulative_distribution_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x1, n1, theta, ...) pnorm(x1, mean = sqrt(n1) * theta, sd = 1) )

setMethod("quantile", signature("DataDistribution"),
          function(x, probs, n1, theta, ...) qnorm(probs, mean = sqrt(n1) * theta, sd = 1) )
