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
#' The logical option \code{two_armed} allows to decide whether a one-arm or
#' a two-arm (the default) design should be computed. In the case of a two-arm
#' design the function \code{\link{plot}} and \code{\link{summary}} describe
#' the sample size per group.
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
#' The option \code{two_armed} can be set to decide whether a one-arm or a
#' two-arm design should be computed.
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
          function(dist, x, n, theta, ...)
              apply(cbind(x, n, theta), 1, function(y)
                  stats::dnorm(y[1], mean = sqrt(y[2]) *
                                   ifelse(dist@two_armed == T, y[3] / sqrt(2), y[3]),
                               sd = 1)
              )
          )


#' @rdname NormalDataDistribution-class
#' @export
setMethod("cumulative_distribution_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...)
              apply(cbind(x, n, theta), 1, function(y)
                  stats::pnorm(y[1], mean = sqrt(y[2]) *
                                   ifelse(dist@two_armed == T, y[3] / sqrt(2), y[3]),
                               sd = 1)
              )
          )



#' @param probs vector of probabilities
#' @rdname NormalDataDistribution-class
#' @export
setMethod("quantile", signature("Normal"),
          function(dist, x, probs, n, theta, ...)
              apply(cbind(probs, n, theta), 1, function(y)
                  stats::qnorm(y[1], mean = sqrt(y[2]) *
                                   ifelse(dist@two_armed == T, y[3] / sqrt(2), y[3]),
                               sd = 1)
              )
          )
