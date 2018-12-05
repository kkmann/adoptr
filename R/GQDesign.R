#' Gaussian-quadrature based implementation of \code{Design}
#'
#' Implements a class for two-stage designs where n2 and c2 are approximated
#' by linear interpolation of pivot points and integration is handled via
#' Gaussian quadrature rules.
#' The pivots are chosen to be the evaluation points of the selected
#' Gaussian quadrature rule and guarantee a stable evaluation during
#' optimization.
#'
#' @slot rule data.frame, the integration rule, result of a call to
#'     \code{gaussquad::legendre.quadrature.rules(order)[[order]])} during
#'     construction
#'
#' @exportClass GQDesign
setClass("GQDesign", representation(
        n1        = "numeric",
        c1f       = "numeric",
        c1e       = "numeric",
        n2_pivots = "numeric",
        c2_pivots = "numeric",
        rule      = "data.frame"
    ),
    contains = "Design")

#' Constructor
#'
#' [TODO: implement a constructor taking in n2 and c2 functions?]
#'
#' @param n1 stage-one sample size
#' @param c1f early stopping for futility boundary
#' @param c1e early stopping for efficacy boundary
#' @param n2_pivots vector of length order giving the values of n2 at the
#'     pivot points of the Gaussian quadrature rule [TODO: these are not available during construction]
#' @param c2_pivots vector of length order giving the values of c2 at the
#'     pivot points of the Gaussian quadrature rule [TODO: these are not available during construction]
#' @param order order (i.e. number of pivot points in the interior of [c1f, c1e])
#'     of the Gaussian quadrature rule to use for integration
#'
#' @rdname GQDesign-class
#' @export
GQDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots, order) {
    if (length(n2_pivots) != order | length(c2_pivots) != order )
        stop("length of pivot vectors does not fit")
    new("GQDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
        c2_pivots = c2_pivots,
        rule = gaussquad::legendre.quadrature.rules(order)[[order]])
}

#' @param object GQDesign to update with given parameter values.
#'
#' @describeIn GQDesign return a new design object with parameters given by a
#'     single numeric vector of length 2*order + 3, inverse of \code{as.numeric(design)}.
setMethod("update", signature("GQDesign"),
    function(object, params, ...) {
        if( ((length(params) - 3) / 2) != nrow(object@rule))
            stop("parameter length does not fit")
        new("GQDesign",
            n1  = params[1],
            c1f = params[2],
            c1e = params[3],
            n2_pivots = params[4:(3 + nrow(object@rule))],
            c2_pivots = params[(4 + nrow(object@rule)):(length(params))],
            rule = object@rule)
    })


#' @param d object of class GQDesign
#' @describeIn GQDesign stage-one sample size
setMethod("n1", signature("GQDesign"), function(d, ...) d@n1)


#' @describeIn GQDesign get the pivots points (knots) of the Gaussian quadrature
#'     rule.
setMethod("get_knots", signature("GQDesign"),
    function(d, ...){
        h <- (d@c1e - d@c1f) / 2
        legendre_knots <- h * d@rule$x + (h + d@c1f)
        # need to reorder pivots (sored in reverse order ...)
        return(legendre_knots[nrow(d@rule):1])
    })


#' @param x1 stage-one outcome
#' @describeIn GQDesign stage-two sample size by linear interpolation of n2_pivots
setMethod("n2", signature("GQDesign", "numeric"),
    function(d, x1, ...) ifelse(x1 < d@c1f | x1 > d@c1e, 0, 1) *
        pmax(0, stats::approx(get_knots(d), d@n2_pivots, xout = x1, method = "linear", rule = 2)$y) )


#' @describeIn GQDesign stage-two critical value by linear interpolation of n2_pivots
setMethod("c2", signature("GQDesign", "numeric"),
    function(d, x1, ...) stats::approx(get_knots(d), d@c2_pivots, xout = x1, method = "linear", rule = 2)$y *
        ifelse(x1 < d@c1f, Inf, 1) * ifelse(x1 > d@c1e, -Inf, 1) )


#' @param x GQDesign
#' @describeIn GQDesign convert GQDesign to vector of tunable parameters via
#'     \code{c(x@n1, x@c1f, x@c1e, x@n2_pivots, x@c2_pivots)}; inverse of update();
setMethod("as.numeric", signature("GQDesign"),
        function(x, ...) c(x@n1, x@c1f, x@c1e, x@n2_pivots, x@c2_pivots))

# not user facing!
setMethod(".evaluate", signature("IntegralScore", "GQDesign"),
    function(s, design, ...) {
        # use design specific implementation tailored to this particular
        # implementation (Gauss Quadrature N points here)
        poef <- predictive_cdf(s@cs@distribution, s@cs@prior, design@c1f, n1(design))
        poee <- 1 - predictive_cdf(s@cs@distribution, s@cs@prior, design@c1e, n1(design))
        # continuation region
        integrand   <- function(z1) evaluate(s@cs, design, z1, ...) *
            predictive_pdf(s@cs@distribution, s@cs@prior, z1, n1(design), ...)
        mid_section <- gaussquad::legendre.quadrature(
            integrand, design@rule, design@c1f, design@c1e
        )
        # compose
        res <- poef * evaluate( # score is constant on early stopping region (TODO: relax later!)
                s@cs, design,
                design@c1f - sqrt(.Machine$double.eps) # slightly smaller than stopping for futility
            ) +
            mid_section +
            poee * evaluate(
                s@cs, design,
                design@c1e + sqrt(.Machine$double.eps)
            )
        return(res)
    })


# not user facing! we need to redo this whole SMoothness stuff...
setMethod(".evaluate", signature("Smoothness_n2", "GQDesign"),
          function(s, design, ...) mean((diff(design@n2_pivots) / diff(get_knots(design)))^2) )
