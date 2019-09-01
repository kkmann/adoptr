#' Univariate prior on model parameter
#'
#' A \code{Prior} object represents a prior distribution on the single model
#' parameter of a \code{\link{DataDistribution}} class
#' object.
#' Together a prior and data-distribution specify the class of the joint
#' distribution of the test statisic, X, and its parameter, theta.
#' Currently, \pkg{adoptr} only allows simple models with a single parameter.
#' Implementations for \link{PointMassPrior} and \link{ContinuousPrior} are available.
#'
#' @examples
#' disc_prior <- PointMassPrior(c(0.1, 0.25), c(0.4, 0.6))
#'
#' cont_prior <- ContinuousPrior(
#'   pdf     = function(x) dnorm(x, mean = 0.3, sd = 0.2),
#'   support = c(-2, 3)
#' )
#'
#'
#' @seealso For the available methods, see \code{\link{bounds}},
#'   \code{\link{expectation}}, \code{\link{condition}}, \code{\link{predictive_pdf}},
#'   \code{\link{predictive_cdf}}, \code{\link{posterior}}
#'
#' @details For an example on working with priors, see
#'    \href{https://kkmann.github.io/adoptr/articles/working-with-priors.html}{here}.
#'
#' @aliases Prior
#' @exportClass Prior
setClass("Prior", representation(label = "character"))




setMethod("show", signature(object = "Prior"), function(object) {
    cat(print(object), "\n")
})





#' Get support of a prior or data distribution
#'
#' \code{bounds()} returns the range of the support of a prior or data distribution.
#'
#' @template dist
#' @template dotdotdot
#'
#' @return \code{numeric} of length two, \code{c(lower, upper)}
#'
#' @export
setGeneric("bounds", function(dist, ...) standardGeneric("bounds"))


#' Expected value of a function
#'
#' Computes the expected value of a vectorized, univariate function \code{f}
#' with respect to a distribution \code{dist}.
#' I.e., \ifelse{html}{\out{E[f(X)]}}{\eqn{\boldsymbol{E}\big[f(X)\big]}{E[f(X)]}}.
#'
#' @param    f     a univariate function, must be vectorized
#' @template dist
#' @template dotdotdot
#'
#' @return \code{numeric}, expected value of \code{f} with respect to \code{dist}
#'
#' @export
setGeneric("expectation", function(dist, f, ...) standardGeneric("expectation"))


#' Condition a prior on an interval
#'
#' Restrict an object of class \code{\link{Prior}} to a sub-interval and
#' re-normalize the PDF.
#'
#' @template dist
#' @param interval length-two numeric vector giving the parameter interval to
#'     condition on
#' @template dotdotdot
#'
#' @return conditional \code{\link{Prior}} on given interval
#'
#' @export
setGeneric("condition", function(dist, interval, ...) standardGeneric("condition"))


#' Predictive PDF
#'
#' \code{predictive_pdf()} evaluates the predictive PDF of the model specified
#' by a \code{\link{DataDistribution}} \code{dist} and
#' \code{\link{Prior}} at the given stage-one outcome.
#'
#' @template dist
#' @template prior
#' @template x1
#' @template n1
#' @template dotdotdot
#'
#' @return \code{numeric}, value of the predictive PDF
#'
#' @export
setGeneric("predictive_pdf", function(dist, prior, x1, n1, ...) standardGeneric("predictive_pdf"))


#' Predictive CDF
#'
#' \code{predictive_cdf()} evaluates the predictive CDF of the model specified
#' by a \code{\link{DataDistribution}} \code{dist} and
#' \code{\link{Prior}} at the given stage-one outcome.
#'
#' @template dist
#' @template prior
#' @template x1
#' @template n1
#' @template dotdotdot
#'
#' @return \code{numeric}, value of the predictive CDF
#'
#' @export
setGeneric("predictive_cdf", function(dist, prior, x1, n1, ...) standardGeneric("predictive_cdf"))


#' Compute posterior distribution
#'
#' Return posterior distribution given observing stage-one outcome.
#'
#' @template dist
#' @template prior
#' @template x1
#' @template n1
#' @template dotdotdot
#'
#' @return Object of class \code{\link{Prior}}
#'
#' @export
setGeneric("posterior", function(dist, prior, x1, n1, ...) standardGeneric("posterior"))
