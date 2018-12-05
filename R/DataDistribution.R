setClass("DataDistribution")

setGeneric("probability_density_function", function(dist, x, n, theta, ...) standardGeneric("probability_density_function"))

setGeneric("cumulative_distribution_function", function(dist, x, n, theta, ...) standardGeneric("cumulative_distribution_function"))

setMethod("quantile", signature("DataDistribution"), function(x, probs, n, theta, ...) stop("not implemented"))







#' Normal data distribution
#'
#' Implements a normal data distribution for z-values given an observed z-value
#' and stage size.
#' Standard deviation is 1 and mean \eqn{\theta\sqrt n} where
#' \eqn{\theta} is the standardized effect size.
#'
#' @exportClass Normal
setClass("Normal", representation(
    dummy = "logical" # needed to make this non-abstract
    ),
    contains = "DataDistribution")


#' Constructor
#'
#' @rdname Normal-class
#' @export
Normal <- function() new("Normal", dummy = FALSE)


#' @param x observed z value
#' @param n stage size
#' @describeIn Normal PDF
#' @export
setMethod("probability_density_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) stats::dnorm(x, mean = sqrt(n) * theta, sd = 1) )


#' @describeIn Normal CDF
#' @export
setMethod("cumulative_distribution_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) stats::pnorm(x, mean = sqrt(n) * theta, sd = 1) )


#' @param probs vector of probablities
#' @describeIn Normal quantile function
#' @export
setMethod("quantile", signature("DataDistribution"),
          function(x, probs, n, theta, ...) stats::qnorm(probs, mean = sqrt(n) * theta, sd = 1) )
