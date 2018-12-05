#' Point-mass prior
#'
#' \code{PointMassPrior} is a generic class for representing a univariate prior
#' over a discrete set of points with positive probability mass.
#' All methods for prior distributions [TODO link] are specialized to
#' the discrete case.
#'
#' @slot theta numeric vector of pivot points (parameter values with positive
#'     prior mass)
#' @slot mass numeric vector of the same length as theta - corresponding
#'     probability masses (must sum to one!)
#'
#' @exportClass PointMassPrior
setClass("PointMassPrior", representation(
        theta = "numeric",
        mass  = "numeric"
    ),
    contains = "Prior")



#' @param theta cf. slot 'theta'
#' @param mass cf. slot 'mass'
#'
#' @rdname PointMassPrior-class
#' @export
PointMassPrior <- function(theta, mass) {
    if (sum(mass) != 1)
        stop("mass must sum to one")
    new("PointMassPrior", theta = theta, mass = mass)
}



#' Get support of prior
#'
#' \code{bounds(dist, ...)} simply returns \code{range(dist@theta)}
#'
#' @rdname PointMassPrior-class
#' @export
setMethod("bounds", signature("PointMassPrior"),
    function(dist, ...) range(dist@theta))


#' @rdname PointMassPrior-class
#' @export
setMethod("expectation", signature("PointMassPrior", "function"),
    function(dist, f, ...) sum(dist@mass * sapply(dist@theta, f, ...)) )


#' @rdname PointMassPrior-class
#' @export
setMethod("condition", signature("PointMassPrior", "numeric"),
    function(dist, interval, ...) {
        if (length(interval) != 2)
            stop("interval must be of length 2")
        if (any(!is.finite(interval)))
            stop("interval must be finite")
        if (diff(interval) < 0)
            stop("interval[2] must be larger or equal to interval[1]")
        epsilon <- sqrt(.Machine$double.eps)
        # find indices of pivots wihtin interval (up to machine precision!)
        idx <- (interval[1] - dist@theta <= epsilon) & (dist@theta - interval[2] <= epsilon)
        # re-normalize and return
        return(PointMassPrior(
            dist@theta[idx], dist@mass[idx] / sum(dist@mass[idx])
        ))
    })


#' @rdname PointMassPrior-class
#' @export
setMethod("predictive_pdf", signature("DataDistribution", "PointMassPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        k   <- length(prior@theta)
        res <- numeric(length(x1))
        for (i in 1:k) {
            res <- res + prior@mass[i] * probability_density_function(dist, x1, n1, prior@theta[i]) # must be implemented
        }
        return(res)
    })


#' @rdname PointMassPrior-class
#' @export
setMethod("predictive_cdf", signature("DataDistribution", "PointMassPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        k   <- length(prior@theta)
        res <- numeric(length(x1))
        for (i in 1:k) {
            res <- res + prior@mass[i] * cumulative_distribution_function(dist, x1, n1, prior@theta[i])
        }
        return(res)
    })


#' @rdname PointMassPrior-class
#' @export
setMethod("posterior", signature("DataDistribution", "PointMassPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        mass <- prior@mass * sapply(prior@theta, function(theta) probability_density_function(dist, x1, n1, theta))
        mass <- mass / sum(mass) # normalize
        return(PointMassPrior(prior@theta, mass))
    })
