#' Abstract univariate prior on model parameter
#'
#' \code{Prior} is an abstract class specifying the interface to be implemented
#' by an prior class.
#' Implementations for \link{PointMassPrior} and \link{ContinuousPrior} are available.
#'
#' @exportClass Prior
setClass("Prior")


#' Get support of prior
#'
#' \code{bounds()} returns the range of the support of a prior [TODO: rename support()?]
#'
#' @param dist a distribution object of class \code{Prior} or \code{DataDistribution}
#' @param ... further optional parameters
#'
#' @rdname Prior-class
#' @export
setGeneric("bounds", function(dist, ...) standardGeneric("bounds"))


#' Expected value of a function
#'
#' \code{expectation()} evaluates the expected value of a function \code{f}
#'     depending on only the parameter value with respect to the prior
#'     distribution \code{dist}.
#'
#' @param f a function of the single model parameter, must be vectorized, i.e.
#'     accept a vector of paramter values and return the vector of function values.
#'
#' @rdname Prior-class
#' @export
setGeneric("expectation", function(dist, f, ...) standardGeneric("expectation"))


#' Condition prior on interval
#'
#' \code{condition()} returns a new object of class prior resulting from
#'     conditioning on the given interval.
#'
#' @param interval length two numeric vector giving the parameter interval to
#'     condition on
#'
#' @rdname Prior-class
#' @export
setGeneric("condition", function(dist, interval, ...) standardGeneric("condition"))


#' Predictive PDF
#'
#' \code{predictive_pdf()} evaluates the predictive pdf of the model specified
#'     by a data distribution dist and prior at he given stage-one outcome
#'     (returns PMF for discrete priors).
#'
#' @param prior a distribution object of class \code{Prior}
#' @param x1 stage-one outcome
#' @param n1 stage-one sample size
#'
#' @rdname Prior-class
#' @export
setGeneric("predictive_pdf", function(dist, prior, x1, n1, ...) standardGeneric("predictive_pdf"))


#' Predictive CDF
#'
#' \code{predictive_cdf()} evaluates the predictive cdf of the model specified
#'     by a data distribution dist and prior at he given stage-one outcome.
#'
#' @rdname Prior-class
#' @export
setGeneric("predictive_cdf", function(dist, prior, x1, n1, ...) standardGeneric("predictive_cdf"))


#' Compute the posterior distribution
#'
#' \code{posterior()} returns a new prior distribution [TODO: okay, maybe rethink naming?]
#'     resulting from conditioning on a given stage-one outcome.
#'
#' @rdname Prior-class
#' @export
setGeneric("posterior", function(dist, prior, x1, n1, ...) standardGeneric("posterior"))
