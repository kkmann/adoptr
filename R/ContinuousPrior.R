setClass("ContinuousPrior", representation(
        pdf     = "function", # actual prior pdf
        support = "numeric"   # the compact support of the prior
    ),
    contains = "Prior")

ContinuousPrior <- function(pdf, support) {
    if (length(support) != 2)
        stop("support must be of length 2")
    if (any(!is.finite(support)))
        stop("support must be finite")
    if (diff(support) <= 0)
        stop("support[2] must be larger (not equal) to support[1]")
    if (abs(stats::integrate(pdf, support[1], support[2], abs.tol = .0005)$value - 1) > .001)
        stop("pdf must integrate to one!")
    new("ContinuousPrior", pdf = pdf, support = support)
}

setMethod("bounds", signature("ContinuousPrior"),
    function(dist, ...) dist@support)

setMethod("expectation", signature("ContinuousPrior", "function"),
    function(dist, f, rel.tol = .001, ...) {
        stats::integrate(
            function(theta) f(theta) * dist@pdf(theta),
            dist@support[1], dist@support[2], rel.tol = rel.tol
        )$value
    })

setMethod("condition", signature("ContinuousPrior", "numeric"),
    function(dist, interval, ...) {
        if (length(interval) != 2)
            stop("interval must be of length 2")
        if (any(!is.finite(interval)))
            stop("interval must be finite")
        if (diff(interval) < 0)
            stop("interval[2] must be larger or equal to interval[1]")
        # compute new normalizing constants
        z <- stats::integrate(dist@pdf, interval[1], interval[2], abs.tol = .00001)$value
        ContinuousPrior(
            function(theta) dist@pdf(theta)/z,
            interval
        )
    })

setMethod("predictive_pdf", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, k = 33, ...) {
        # TODO: use Gaussian Quadrature for pivots!
        piv <- seq(prior@support[1], prior@support[2], length.out = k)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        res <- numeric(length(x1))
        for (i in 1:k) {
            res <- res + mass[i] * probability_density_function(dist, x1, n1, piv[i])
        }
        return(res)
    })

setMethod("predictive_cdf", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, k = 33, ...) {
        # TODO: use Gaussian Quadrature for pivots!
        piv <- seq(prior@support[1], prior@support[2], length.out = k)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        res <- numeric(length(x1))
        for (i in 1:k) {
            res <- res + mass[i] * cumulative_distribution_function(dist, x1, n1, piv[i])
        }
        return(res)
    })

setMethod("posterior", signature("DataDistribution", "ContinuousPrior", "numeric"),
    function(dist, prior, x1, n1, ...) { # careful, assumed null hypothesis = 0
        if (length(x1) != 1)
            stop("no vectorized version in x1")
        # TODO: use Gaussian Quadrature
        prop_pdf <- function(theta) probability_density_function(dist, x1, n1, theta) * prior@pdf(theta)
        z <- stats::integrate(prop_pdf, prior@support[1], prior@support[2], abs.tol = .00001)$value
        ContinuousPrior(
            function(theta) prop_pdf(theta)/z, prior@support
        )
    })
