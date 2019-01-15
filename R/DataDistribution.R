#' Data distributions
#'
#' \code{DataDistribution} is an abstract class used to represent the distribution
#' of a sufficient statistic \code{x} given a sample size \code{n} and a
#' single parameter value \code{theta}.
#'
#' This abstraction layer allows representation of t-distributions
#' (unknown variance), normal distribution (known variance), and normal
#' approximation of a binary endpoint.
#' Currently, only the normal case is implemented with \code{\link{Normal-class}}.
#'
#' @slot two_armed Should a two-armed design be used?
#'
#' @template DataDistributionTemplate
#'
#' @exportClass DataDistribution
setClass("DataDistribution", representation(
    two_armed = "logical")
)


#' @param probs numeric vector of probabilities
#' @describeIn DataDistribution quantile function of the respective distribution.
#' @export
setMethod("quantile", signature("DataDistribution"), function(dist, x, probs, n, theta, ...) stop("not implemented"))


#' @rdname DataDistribution-class
#' @export
setGeneric("probability_density_function", function(dist, x, n, theta, ...) standardGeneric("probability_density_function"))

#' @describeIn DataDistribution probability density function given outcome,
#'     sample size and parameter; must be implemented.
#' @export
setMethod("probability_density_function", signature("DataDistribution", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) stop("not implemented"))


#' @rdname DataDistribution-class
#' @export
setGeneric("cumulative_distribution_function", function(dist, x, n, theta, ...) standardGeneric("cumulative_distribution_function"))

#' @describeIn DataDistribution cumulative distribution function given outcome,
#'     sample size and parameter; must be implemented.
#' @export
setMethod("cumulative_distribution_function", signature("DataDistribution", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) stop("not implemented"))





#' Normal data distribution
#'
#' Implements a normal data distribution for z-values given an observed z-value
#' and stage size.
#' Standard deviation is 1 and mean \eqn{\theta\sqrt n} where
#' \eqn{\theta} is the standardized effect size.
#' See \code{\link{DataDistribution-class}} for more details.
#'
#' @template DataDistributionTemplate
#'
#' @rdname NormalDataDistribution-class
#' @exportClass Normal
setClass("Normal", representation(
    two_armed = "logical"
    ),
    contains = "DataDistribution")

#' @param two_armed s. slot
#'
#' @rdname NormalDataDistribution-class
#' @export
Normal <- function(two_armed = TRUE) new("Normal", two_armed = two_armed)


#' @rdname NormalDataDistribution-class
#' @export
setMethod("probability_density_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...){
              theta <- ifelse(dist@two_armed == T, theta / sqrt(2), theta)
              return(stats::dnorm(x, mean = sqrt(n) * theta, sd = 1))
              }
          )


#' @rdname NormalDataDistribution-class
#' @export
setMethod("cumulative_distribution_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...){
              theta <- ifelse(dist@two_armed == T, theta / sqrt(2), theta)
              return(stats::pnorm(x, mean = sqrt(n) * theta, sd = 1))
              }
          )



#' @param probs vector of probabilities
#' @rdname NormalDataDistribution-class
#' @export
setMethod("quantile", signature("Normal"),
          function(dist, x, probs, n, theta, ...){
              theta <- ifelse(dist@two_armed == T, theta / sqrt(2), theta)
              return(stats::qnorm(probs, mean = sqrt(n) * theta, sd = 1))
              }
          )
