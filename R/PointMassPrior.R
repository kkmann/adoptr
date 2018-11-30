setClass("PointMassPrior", representation(
        theta = "numeric",
        mass  = "numeric"
    ),
    contains = "Prior")

PointMassPrior <- function(theta, mass) {
    if (sum(mass) != 1)
        stop("mass must sum to one")
    new("PointMassPrior", theta = theta, mass = mass)
}

setMethod("bounds", signature("PointMassPrior"),
    function(prior, ...) range(prior@theta))

setMethod("expectation", signature("PointMassPrior", "function"),
    function(prior, f, ...) sum(prior@mass * sapply(prior@theta, f, ...)) )

setMethod("condition", signature("PointMassPrior", "numeric"),
    function(prior, interval, ...) {
        if (length(interval) != 2)
            stop("interval must be of length 2")
        if (any(!is.finite(interval)))
            stop("interval must be finite")
        if (diff(interval) < 0)
            stop("interval[2] must be larger or equal to interval[1]")
        epsilon <- sqrt(.Machine$double.eps)
        # find indices of pivots wihtin interval (up to machine precision!)
        idx <- (interval[1] - prior@theta <= epsilon) & (prior@theta - interval[2] <= epsilon)
        # re-normalize and return
        return(PointMassPrior(
            prior@theta[idx], prior@mass[idx] / sum(prior@mass[idx])
        ))
    })

setMethod("predictive_pdf", signature("PointMassPrior", "numeric"),
    function(prior, z1, n1, ...) {
        k   <- length(prior@theta)
        res <- numeric(length(z1))
        for (i in 1:k) { # TODO: careful, assumed null hypothesis = 0
            res <- res + prior@mass[i] * dnorm(z1, mean = sqrt(n1) * prior@theta[i], sd = 1)
        }
        return(res)
    })

setMethod("predictive_cdf", signature("PointMassPrior", "numeric"),
    function(prior, z1, n1, ...) {
        k   <- length(prior@theta)
        res <- numeric(length(z1))
        for (i in 1:k) { # TODO: careful, assumed null hypothesis = 0
            res <- res + prior@mass[i] * pnorm(z1, mean = sqrt(n1) * prior@theta[i], sd = 1)
        }
        return(res)
    })

setMethod("posterior", signature("PointMassPrior", "numeric"),
    function(prior, z1, n1, ...) { # careful, assumed null hypothesis = 0
        mass <- prior@mass * dnorm(z1, mean = sqrt(n1) * prior@theta, sd = 1)
        mass <- mass / sum(mass) # normalize
        return(PointMassPrior(prior@theta, mass))
    })
