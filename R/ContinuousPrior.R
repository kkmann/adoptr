#' Continuous univariate prior distributions
#'
#' \code{ContinuousPrior} is a sub-class of \code{\link{Prior}} implementing
#' a generic representation of continuous prior distributions over a compact
#' interval on the real line.
#'
#' @slot pdf     cf. parameter 'pdf'
#' @slot support cf. parameter 'support'
#' @slot pivots normalized pivots for integration rule (in [-1, 1])
#'     the actual pivots are scaled to the support of the prior
#' @slot weights weights of of integration rule at \code{pivots} for
#'     approximating integrals over \code{delta}
#'
#' @seealso Discrete priors are supported via \code{\link{PointMassPrior}}
#'
#' @aliases ContinuousPrior
#' @exportClass ContinuousPrior
setClass("ContinuousPrior", representation(
        pdf     = "function",
        support = "numeric",
        pivots  = "numeric",
        weights = "numeric"
        ),
    contains = "Prior")


#' @param pdf                 vectorized univariate PDF function
#' @param support             numeric vector of length two with the bounds of
#'     the compact interval on which the pdf is positive.
#' @template order
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
                            order = 10,
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

    rule   <- GaussLegendreRule(order)
    h      <- (support[2] - support[1]) / 2
    pivots <- h * rule$nodes + (h + support[1])


    new("ContinuousPrior", pdf = pdf, support = support,
        pivots = pivots, weights = rule$weights, label = label)
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
#' @rdname expectation
#' @export
setMethod("expectation", signature("ContinuousPrior", "function"),
    function(dist, f, ...) {
        gauss_quad(f(dist@pivots) * dist@pdf(dist@pivots), dist@support[1], dist@support[2], dist@weights)
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
        if (length(interval) != 2) stop("interval must be of length 2")
        if (any(!is.finite(interval))) stop("interval must be finite")
        if (diff(interval) < 0) stop("interval[2] must be larger or equal to interval[1]")
        interval[1] <- max(interval[1], dist@support[1])
        interval[2] <- min(interval[2], dist@support[2])
        if (diff(interval) < 0) stop("resulting interval is empty")
        z       <- stats::integrate(dist@pdf, interval[1], interval[2], abs.tol = .001)$value
        new_pdf <- function(theta) ifelse(interval[1] <= theta & theta <= interval[2], dist@pdf(theta)/z, 0)
        h_old   <- (dist@support[2] - dist@support[1]) / 2
        nodes   <- (dist@pivots - (h_old + dist@support[1])) / h_old
        h_new   <- (interval[2] - interval[1]) / 2
        pivots  <- h_new * nodes + (h_new + interval[1])

        new("ContinuousPrior", pdf = new_pdf, support = interval,
            pivots = pivots, weights = dist@weights, label = dist@label)
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
    function(dist, prior, x1, n1, k = 10*(prior@support[2] - prior@support[1]) + 1, ...) {
        piv  <- seq(prior@support[1], prior@support[2], length.out = k)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        grid <- expand.grid(x1 = x1, piv = piv)
        (matrix(probability_density_function(dist, grid$x1, n1, grid$piv), nrow = length(x1)) %*% mass)[, 1]
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
    function(dist, prior, x1, n1, k = 10*(prior@support[2] - prior@support[1]) + 1, ...) {
        piv  <- seq(prior@support[1], prior@support[2], length.out = k)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        grid <- expand.grid(x1 = x1, piv = piv)
        (matrix(cumulative_distribution_function(dist, grid$x1, n1, grid$piv), nrow = length(x1)) %*% mass)[, 1]
    })


#' @examples
#' tmp <- ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4))
#' posterior(Normal(), tmp, 2, 20)
#'
#' @rdname posterior
#' @export
setMethod("posterior", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        if (length(x1) != 1) stop("no vectorized version in x1")
        prop_pdf <- function(theta) probability_density_function(dist, x1, n1, theta) * prior@pdf(theta)
        z        <- gauss_quad(prop_pdf(prior@pivots), prior@support[1], prior@support[2], prior@weights)
        new("ContinuousPrior",
            pdf     = function(theta) prop_pdf(theta) / z,
            support = prior@support,
            pivots  = prior@pivots,
            weights = prior@weights,
            label   = prior@label
        )
    })



setMethod("print", signature('ContinuousPrior'), function(x, ...) {
    name <- if (!is.na(x@label)) x@label else class(x)[1]
    glue::glue("{name}<[{x@support[1]},{x@support[2]}]>")
})
