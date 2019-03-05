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
#' design all sample sizes are per group.
#'
#' @slot two_armed Should a two-armed design be used?
#'
#' @template DataDistributionTemplate
#'
#' @exportClass DataDistribution
setClass("DataDistribution", representation(
    two_armed = "logical")
)



#' @rdname DataDistribution-class
#' @export
setGeneric("probability_density_function", function(dist, x, n, theta, ...) standardGeneric("probability_density_function"))



#' @rdname DataDistribution-class
#' @export
setGeneric("cumulative_distribution_function", function(dist, x, n, theta, ...) standardGeneric("cumulative_distribution_function"))





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
          function(dist, x, n, theta, ...) {
              if (dist@two_armed) {
                  theta <- theta / sqrt(2)
              }
              stats::dnorm(x, mean = sqrt(n) * theta, sd = 1)
          })



#' @rdname NormalDataDistribution-class
#' @export
setMethod("cumulative_distribution_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              if (dist@two_armed) {
                  theta <- theta / sqrt(2)
              }
              stats::pnorm(x, mean = sqrt(n) * theta, sd = 1)
          })



#' @param probs vector of probabilities
#' @rdname NormalDataDistribution-class
#' @export
setMethod("quantile", signature("Normal"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
              if (x@two_armed) {
                  theta <- theta / sqrt(2)
              }
              stats::qnorm(probs, mean = sqrt(n) * theta, sd = 1)
          })



#' @rdname NormalDataDistribution-class
#'
#' @param object design to simulate from
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("Normal", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
              fct <- 1
              if (object@two_armed)
                  fct <- 1 / sqrt(2)

              if (!is.null(seed))
                  set.seed(seed)

              stats::rnorm(nsim, mean = fct * sqrt(n) * theta, sd = 1)
          })


#' @rdname NormalDataDistribution-class
#'
#' @param object object of class \code{Normal}
#' @export
setMethod("show", signature(object = "Normal"),
          function(object) cat(class(object)[1]))

