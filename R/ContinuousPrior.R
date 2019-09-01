#' Continuous univariate prior distributions
#'
#' \code{ContinuousPrior} is a sub-class of \code{\link{Prior}} implementing
#' a generic representation of continuous prior distributions over a compact
#' interval on the real line.
#'
#' @slot pdf     cf. parameter 'pdf'
#' @slot support cf. parameter 'support'
#'
#' @seealso Discrete priors are supported via \code{\link{PointMassPrior}}
#'
#' @aliases ContinuousPrior
#' @exportClass ContinuousPrior
setClass("ContinuousPrior", representation(
        pdf     = "function",
        support = "numeric"
        ),
    contains = "Prior")


#' @param pdf                 vectorized univariate PDF function
#' @param support             numeric vector of length two with the bounds of
#'     the compact interval on which the pdf is positive.
#' @template tighten_support
#' @template check_normalization
#' @template label
#'
#' @examples
#' ContinuousPrior(function(x) 2*x, c(0, 1))
#'
#' @rdname ContinuousPrior-class
#' @export
ContinuousPrior <- function(pdf,
                            support,
                            label = NA_character_,
                            tighten_support = FALSE,
                            check_normalization = TRUE) {
    if (length(support) != 2)
        stop("support must be of length 2")

    if (any(!is.finite(support)))
        stop("support must be finite")

    if (diff(support) <= 0)
        stop("support[2] must be larger (not equal) to support[1]")

    if (check_normalization) {
        if (abs(stats::integrate(pdf, support[1], support[2], abs.tol = .00001)$value - 1) > .001)
            stop("pdf must integrate to one!")
    }

    if (tighten_support) {
        eps <- .001
        while (pdf(support[1] + eps) < .Machine$double.eps^2) {
            support[1] <- support[1] + eps
        }
        while (pdf(support[2] - eps) < .Machine$double.eps^2) {
            support[2] <- support[2] - eps
        }
        norm    <- stats::integrate(pdf, support[1], support[2], abs.tol = .00001)$value
        pdf_old <- pdf
        pdf     <- function(theta) pdf_old(theta) / norm
    }

    new("ContinuousPrior", pdf = pdf, support = support, label = label)
}


#' @examples
#' bounds(ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4)))
#' # > 0.2 0.4
#'
#' @rdname bounds
#' @export
setMethod("bounds", signature("ContinuousPrior"),
    function(dist, ...) dist@support)


#' @examples
#' expectation(
#'     ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4)),
#'     identity
#' )
#' # > 0.3
#'
#' @param rel.tol \code{numeric}, relative tolerance used in adaptive gaussian
#'     quadrature to evaluate the integral
#'
#' @rdname expectation
#' @export
setMethod("expectation", signature("ContinuousPrior", "function"),
    function(dist, f, rel.tol = .001, ...) {
        stats::integrate(
            function(theta) f(theta) * dist@pdf(theta),
            dist@support[1], dist@support[2], rel.tol = rel.tol
        )$value
    })


#' @examples
#' tmp <- condition(
#'     ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4)),
#'     c(.3, .5)
#' )
#' bounds(tmp) # c(.3, .4)
#'
#' @rdname condition
#' @export
setMethod("condition", signature("ContinuousPrior", "numeric"),
    function(dist, interval, ...) {
        if (length(interval) != 2)
            stop("interval must be of length 2")
        if (any(!is.finite(interval)))
            stop("interval must be finite")
        if (diff(interval) < 0)
            stop("interval[2] must be larger or equal to interval[1]")

        interval[1] <- max(interval[1], dist@support[1])
        interval[2] <- min(interval[2], dist@support[2])
        if (diff(interval) < 0)
            stop("resulting interval is empty")
        # compute new normalizing constant
        z <- stats::integrate(dist@pdf, interval[1], interval[2], abs.tol = .00001)$value
        new_pdf <- function(theta) {
            ifelse(interval[1] <= theta & theta <= interval[2],
                dist@pdf(theta) / z,
                0
            )
        }
        ContinuousPrior(
            new_pdf,
            interval
        )
    })


#' @examples
#' tmp <- ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4))
#' predictive_pdf(Normal(), tmp, 2, 20)
#'
#' @param  k number of pivots for crude integral approximation
#'
#' @rdname predictive_pdf
#' @export
setMethod("predictive_pdf", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, k = 33, ...) {
        piv  <- seq(prior@support[1], prior@support[2], length.out = k)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        res  <- numeric(length(x1))
        for (i in 1:k) {
            res <- res + mass[i] * probability_density_function(dist, x1, n1, piv[i])
        }
        return(res)
    })


#' @examples
#' tmp <- ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4))
#' predictive_cdf(Normal(), tmp, 2, 20)
#'
#' @param  k number of pivots for crude integral approximation
#'
#' @rdname predictive_cdf
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


#' @examples
#' tmp <- ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4))
#' posterior(Normal(), tmp, 2, 20)
#'
#' @template tighten_support
#' @template check_normalization
#'
#' @rdname posterior
#' @export
setMethod("posterior", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, tighten_support = FALSE, check_normalization = FALSE, ...) {
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
            tighten_support     = tighten_support,
            check_normalization = check_normalization
        )
    })



setMethod("print", signature('ContinuousPrior'), function(x, ...) {
    name <- if (!is.na(x@label)) x@label else class(x)[1]
    glue::glue("{name}<[{x@support[1]},{x@support[2]}]>")
})
