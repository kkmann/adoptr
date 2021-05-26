#' Data distributions
#'
#' \code{DataDistribution} is an abstract class used to represent the distribution
#' of a sufficient statistic \code{x} given a sample size \code{n} and a
#' single parameter value \code{theta}.
#'
#' This abstraction layer allows the representation of t-distributions
#' (unknown variance), normal distribution (known variance), and normal
#' approximation of a binary endpoint.
#' Currently, the two implemented versions are \code{\link{Normal-class}} and
#' \code{\link{Binomial-class}}.
#'
#' The logical option \code{two_armed} allows to decide whether a one-arm or
#' a two-arm (the default) design should be computed. In the case of a two-arm
#' design all sample sizes are per group.
#'
#' @slot two_armed Logical that indicates if a two-arm design is assumed.
#'
#' @examples
#' normaldist   <- Normal(two_armed = FALSE)
#' binomialdist <- Binomial(rate_control = .25, two_armed = TRUE)
#'
#' @template DataDistributionTemplate
#'
#' @aliases DataDistribution
#' @exportClass DataDistribution
setClass("DataDistribution", representation(
    two_armed = "logical")
)


#' Probability density function
#'
#' \code{probability_density_function} evaluates the probability density
#' function of a specific distribution \code{dist} at a point \code{x}.
#'
#' @template dist
#' @template DataDistributionTemplate
#'
#' @export
setGeneric("probability_density_function", function(dist, x, n, theta, ...) standardGeneric("probability_density_function"))


#' Cumulative distribution function
#'
#' \code{cumulative_distribution_function} evaluates the cumulative distribution
#' function of a specific distribution \code{dist} at a point \code{x}.
#'
#' @template dist
#' @template DataDistributionTemplate
#'
#' @export
setGeneric("cumulative_distribution_function", function(dist, x, n, theta, ...) standardGeneric("cumulative_distribution_function"))


setMethod("show", signature(object = "DataDistribution"), function(object) {
    cat(print(object), "\n")
})
