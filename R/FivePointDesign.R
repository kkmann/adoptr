setClass("FivePointDesign", representation(
        n1        = "numeric",
        c1f       = "numeric",
        c1e       = "numeric",
        n2_pivots = "numeric",
        c2_pivots = "numeric"
    ),
    contains = "Design")

FivePointDesign <- function(n1, c1f, c1e, n2_pivots, c2_pivots) {
    if (length(n2_pivots) != 5 | length(c2_pivots) != 5)
        stop("both pivot vectors must have length 5")
    new("FivePointDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
        c2_pivots = c2_pivots)
}

setMethod("update", signature("FivePointDesign"),
    function(object, params, ...) {
        if (length(params) != 13) # n1, c1f, c1e, 2*5 for n2/c2 pivots
            stop("parameter length does not fit")
        new("FivePointDesign", n1 = params[1], c1f = params[2], c1e = params[3],
          n2_pivots = params[4:(3 + 5)], c2_pivots = params[(4 + 5):13])
    })

setMethod("n1", signature("FivePointDesign"), function(d, ...) d@n1)


setGeneric("get_knots", function(d, ...) standardGeneric("get_knots"))
setMethod("get_knots", signature("FivePointDesign"),
    function(d, ...) seq(d@c1f, d@c1e, length.out = 5))

setMethod("n2", signature("FivePointDesign", "numeric"),
    function(d, z1, ...) ifelse(z1 < d@c1f | z1 > d@c1e, 0, 1) *
        pmax(0, approx(get_knots(d), d@n2_pivots, xout = z1, method = "linear", rule = 2)$y) )

setMethod("c2", signature("FivePointDesign", "numeric"),
    function(d, z1, ...) approx(get_knots(d), d@c2_pivots, xout = z1, method = "linear", rule = 2)$y *
        ifelse(z1 < d@c1f, Inf, 1) * ifelse(z1 > d@c1e, -Inf, 1) )

setMethod("as.numeric", signature("FivePointDesign"),
        function(x, ...) c(x@n1, x@c1f, x@c1e, x@n2_pivots, x@c2_pivots))

setMethod(".eval_specific", signature("IntegralScore", "FivePointDesign"),
    function(s, design, ...) {
        # use design specific implementation tailored to this particular
        # implementation (Newton Cotes 5 points here)
        poef <- predictive_cdf(s@conditional_score@prior, design@c1f, n1(design))
        poee <- 1 - predictive_cdf(s@conditional_score@prior, design@c1e, n1(design))
        # continuation region
        integrand   <- function(z1) eval(s@conditional_score, design, z1, ...) *
            predictive_pdf(s@conditional_score@prior, z1, n1(design), ...)
        weights     <- c(7, 32, 12, 32, 7)
        h           <- (design@c1e - design@c1f)/4
        mid_section <- 2/45 * h * sum(weights * integrand(get_knots(design)))
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

setMethod(".eval_specific", signature("Smoothness_n2", "FivePointDesign"),
          function(s, design, ...) mean(diff(design@n2_pivots)^2) )
