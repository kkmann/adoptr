#' Continuous univariate prior distribution
#'
#' \code{ContinuousPrior} is a generic class for representing a univariate prior
#' over a compact interval of the real numbers.
#'
#' For further details cf. \code{\link{Prior-class}}.
#'
#' @slot pdf function the actual probability density function , only called
#'     in the interior of support therefore need not be set to zero outside of
#'     support. Must be vectorized! [TODO: expose via pdf(prior, theta)!]
#' @slot support numeric vector of length two with the bounds of the compact
#'     interval on which the pdf is positive.
#'     The restriction to compact support simplifies numerical integration.
#'
#' @template PriorTemplate
#'
#' @exportClass ContinuousPrior
setClass("ContinuousPrior", representation(
        pdf     = "function", # actual prior pdf
        support = "numeric"   # the compact support of the prior
        ),
    contains = "Prior")


#' @param pdf cf. slot \code{pdf}
#' @param support cf. slot \code{support}
#' @param tighten_support logical indicating if the support should be tightened
#'     automatically.
#'
#' @describeIn ContinuousPrior-class constructor
#' @export
ContinuousPrior <- function(pdf, support, tighten_support = FALSE) {
    if (length(support) != 2)
        stop("support must be of length 2")
    if (any(!is.finite(support)))
        stop("support must be finite")
    if (diff(support) <= 0)
        stop("support[2] must be larger (not equal) to support[1]")
    if(tighten_support) {
        while(pdf(support[1]) < .Machine$double.eps^2) {
            support[1] <- support[1] + .001
            }
        while(pdf(support[2]) < .Machine$double.eps^2) {
            support[2] <- support[2] - .001
            }
    }
    if (abs(stats::integrate(pdf, support[1], support[2], abs.tol = .00001)$value - 1) > .001)
        stop("pdf must integrate to one!")
    new("ContinuousPrior", pdf = pdf, support = support)
}


#' @rdname ContinuousPrior-class
#' @export
setMethod("bounds", signature("ContinuousPrior"),
    function(dist, ...) dist@support)


#' @param rel.tol relative tolerance used in adaptive gaussian quadrature
#'     to evaluate the integral
#' @rdname ContinuousPrior-class
#' @export
setMethod("expectation", signature("ContinuousPrior", "function"),
    function(dist, f, rel.tol = .001, ...) {
        stats::integrate(
            function(theta) f(theta) * dist@pdf(theta),
            dist@support[1], dist@support[2], rel.tol = rel.tol
        )$value
    })


#' @rdname ContinuousPrior-class
#' @export
setMethod("condition", signature("ContinuousPrior", "numeric"),
    function(dist, interval, ...) {
        if (length(interval) != 2)
            stop("interval must be of length 2")
        if (any(!is.finite(interval)))
            stop("interval must be finite")
        if (diff(interval) < 0)
            stop("interval[2] must be larger or equal to interval[1]")
        # compute new normalizing constant
        z <- stats::integrate(dist@pdf, interval[1], interval[2], abs.tol = .00001)$value
        ContinuousPrior(
            function(theta) dist@pdf(theta)/z, interval
        )
    })


#' @param k number of pivots for crude integral approximation [TODO: this needs to be done properly!, cant use integrate since we wnat this vectorized!]
#' @rdname ContinuousPrior-class
#' @export
setMethod("predictive_pdf", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, k = 33, ...) {
        piv <- seq(prior@support[1], prior@support[2], length.out = k)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        res <- numeric(length(x1))
        for (i in 1:k) {
            res <- res + mass[i] * probability_density_function(dist, x1, n1, piv[i])
        }
        return(res)
    })


#' @rdname ContinuousPrior-class
#' @export
setMethod("predictive_cdf", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, k = 33, ...) {
        piv  <- seq(prior@support[1], prior@support[2], length.out = k)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        res  <- numeric(length(x1))
        for (i in 1:k) {
            res <- res + mass[i] * cumulative_distribution_function(dist, x1, n1, piv[i])
        }
        return(res)
    })


#' @rdname ContinuousPrior-class
#' @export
setMethod("posterior", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        if (length(x1) != 1)
            stop("no vectorized version in x1")
        prop_pdf <- function(theta) {
            probability_density_function(dist, x1, n1, theta) * prior@pdf(theta)
        }
        z <- stats::integrate(
            prop_pdf, prior@support[1], prior@support[2], abs.tol = .00001
        )$value
        ContinuousPrior(
            function(theta) prop_pdf(theta) / z,
            prior@support,
            tighten_support = FALSE
        )
    })


#' @rdname ContinuousPrior-class
#'
#' @param object object of class \code{ContinuousPrior}
#' @export
setMethod("show", signature(object = "ContinuousPrior"),
          function(object) cat(class(object)[1]))

