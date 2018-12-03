setClass("GQDesign", representation(
        n1        = "numeric",
        c1f       = "numeric",
        c1e       = "numeric",
        n2_pivots = "numeric",
        c2_pivots = "numeric",
        nb_knots  = "numeric"
    ),
    contains = "Design")

GQDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots, nb_knots) {
    if (length(n2_pivots) != nb_knots | length(c2_pivots) != nb_knots)
        stop("length of pivot vectors does not fit")
    new("GQDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
        c2_pivots = c2_pivots, nb_knots = nb_knots)
}

setMethod("update", signature("GQDesign"),
    function(object, params, ...) {
        #if (length(params) != 13)
        #    stop("parameter length does not fit")
        new("GQDesign",
            n1 = params[1],
            c1f = params[2],
            c1e = params[3],
            n2_pivots = params[4:(3 + object@nb_knots)],
            c2_pivots = params[(4 + object@nb_knots):(length(params))],
            nb_knots = object@nb_knots)
    })

setMethod("n1", signature("GQDesign"), function(d, ...) d@n1)

# Define weights
# Can we define this as global object?
# For runtime issues it is unnecessary to call it in every function call
library(gaussquad)
weights.global <- gaussquad::legendre.quadrature.rules(100)


setGeneric("get_knots", function(d, ...) standardGeneric("get_knots"))
setMethod("get_knots", signature("GQDesign"),
    function(d, ...){
        h <- (d@c1e - d@c1f) / 2
        legendre_knots <- h * weights.global[[d@nb_knots]]$x + (h + d@c1f)
        return(legendre_knots[d@nb_knots:1])
    }
    )

setMethod("n2", signature("GQDesign", "numeric"),
    function(d, z1, ...) ifelse(z1 < d@c1f | z1 > d@c1e, 0, 1) *
        pmax(0, approx(get_knots(d), d@n2_pivots, xout = z1, method = "linear", rule = 2)$y) )

setMethod("c2", signature("GQDesign", "numeric"),
    function(d, z1, ...) approx(get_knots(d), d@c2_pivots, xout = z1, method = "linear", rule = 2)$y *
        ifelse(z1 < d@c1f, Inf, 1) * ifelse(z1 > d@c1e, -Inf, 1) )

setMethod("as.numeric", signature("GQDesign"),
        function(x, ...) c(x@n1, x@c1f, x@c1e, x@n2_pivots, x@c2_pivots))


setMethod(".eval_specific", signature("IntegralScore", "GQDesign"),
    function(s, design, ...) {
        # use design specific implementation tailored to this particular
        # implementation (Gauss Quadrature N points here)
        poef <- predictive_cdf(s@conditional_score@prior, design@c1f, n1(design))
        poee <- 1 - predictive_cdf(s@conditional_score@prior, design@c1e, n1(design))
        # continuation region
        integrand   <- function(z1) eval(s@conditional_score, design, z1, ...) *
            predictive_pdf(s@conditional_score@prior, z1, n1(design), ...)
        weights     <- weights.global[[design@nb_knots]]$w
        h           <- (design@c1e - design@c1f) / 2
        mid_section <- h * sum(weights * integrand(get_knots(design)))
        # compose
        res <- poef * eval( # score is constant on early stopping region (TODO: relax later!)
                s@conditional_score, design,
                design@c1f - sqrt(.Machine$double.eps) # slightly smaller than stopping for futility
            ) +
            mid_section +
            poee * eval(
                s@conditional_score, design,
                design@c1e + sqrt(.Machine$double.eps)
            )
        return(res)
    })

setMethod(".eval_specific", signature("Smoothness_n2", "GQDesign"),
          function(s, design, ...) mean((diff(design@n2_pivots) / diff(get_knots(design)))^2) )
