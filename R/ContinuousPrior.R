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
    if (abs(integrate(pdf, support[1], support[2], abs.tol = .0005)$value - 1) > .001)
        stop("pdf must integrate to one!")
    new("ContinuousPrior", pdf = pdf, support = support)
}

setMethod("bounds", signature("ContinuousPrior"),
    function(prior, ...) prior@support)

setMethod("expectation", signature("ContinuousPrior", "function"),
    function(prior, f, rel.tol = .001, ...) {
        stats::integrate(
            function(theta) f(theta) * prior@pdf(theta),
            prior@support[1], prior@support[2], rel.tol = rel.tol
        )
    })

setMethod("condition", signature("ContinuousPrior", "numeric"),
    function(prior, interval, ...) {
        if (length(interval) != 2)
            stop("interval must be of length 2")
        if (any(!is.finite(interval)))
            stop("interval must be finite")
        if (diff(interval) < 0)
            stop("interval[2] must be larger or equal to interval[1]")
        # compute new normalizing constants
        # TODO: change to non-adaptive quadrature!
        z <- stats::integrate(prior@pdf, interval[1], interval[2], abs.tol = .00001)$value
        ContinuousPrior(
            function(theta) prior@pdf(theta)/z,
            interval
        )
    })

setMethod("predictive_pdf", signature("ContinuousPrior", "numeric"),
    function(prior, z1, n1, k = 33, ...) {
        # TODO: use Gaussian Quadrature for pivots!
        piv <- seq(prior@support[1], prior@support[2], length.out = 100)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        res <- numeric(length(z1))
        for (i in 1:k) { # TODO: careful, assumed null hypothesis = 0
            res <- res + mass[i] * dnorm(z1, mean = sqrt(n1) * piv[i], sd = 1)
        }
        return(res)
    })

setMethod("predictive_cdf", signature("ContinuousPrior", "numeric"),
    function(prior, z1, n1, k = 33, ...) {
        # TODO: use Gaussian Quadrature for pivots!
        piv <- seq(prior@support[1], prior@support[2], length.out = 100)
        mass <- sapply(piv, prior@pdf)
        mass <- mass / sum(mass) # (renormalize!)
        res <- numeric(length(z1))
        for (i in 1:k) { # TODO: careful, assumed null hypothesis = 0
            res <- res + mass[i] * pnorm(z1, mean = sqrt(n1) * piv[i], sd = 1)
        }
        return(res)
    })

setMethod("posterior", signature("ContinuousPrior", "numeric"),
    function(prior, z1, n1, ...) { # careful, assumed null hypothesis = 0
        if (length(z1) != 1)
            stop("no vectorized version in z1")
        # TODO: use Gaussian Quadrature
        prop_pdf <- function(theta) dnorm(z1, mean = sqrt(n1) * theta, sd = 1) * prior@pdf(theta)
        z <- stats::integrate(prop_pdf, interval[1], interval[2], abs.tol = .00001
        )$value
        ContinuousPrior(
            function(theta) prop_pdf(theta)/z, prior@support
        )
    })
