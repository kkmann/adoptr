#' Scores
#'
#' In \code{adoptr} scores are used to assess the performance of a design.
#' This can be done either conditionally on the observed stage-one outcome
#' or unconditionally.
#' Consequently, score objects are either of class \code{ConditionalScore} or
#' \code{UnconditionalScore}.
#'
#' @template s
#' @template design
#' @template dotdotdot
#' @param data_distribution \code{\link{DataDistribution}} object
#' @template prior
#' @template label
#' @template optimization
#' @param subdivisions maximal number of subdivisions when evaluating an integral
#'   score using adaptive quadrature (optimization = FALSE)
#'
#' @details
#' All scores can be evaluated on a design using the \code{evaluate} method.
#' Note that \code{evaluate} requires a third argument \code{x1} for
#' conditional scores (observed stage-one outcome).
#' Any \code{ConditionalScore} can be converted to a \code{UnconditionalScore}
#' by forming its expected value using \code{expected}.
#' The returned unconditional score is of class \code{IntegralScore}.
#'
#' @seealso \code{\link{ConditionalPower}}, \code{\link{ConditionalSampleSize}},
#' \code{\link{composite}}
#'
#' @examples
#' design <- TwoStageDesign(
#'   n1    = 25,
#'   c1f   = 0,
#'   c1e   = 2.5,
#'   n2    = 50,
#'   c2    = 1.96,
#'   order = 7L
#' )
#' prior <- PointMassPrior(.3, 1)
#'
#' # conditional
#' cp <- ConditionalPower(Normal(), prior)
#' expected(cp, Normal(), prior)
#' evaluate(cp, design, x1 = .5)
#'
#' # unconditional
#' power <- Power(Normal(), prior)
#' evaluate(power, design)
#' evaluate(power, design, optimization = TRUE) # use non-adaptive quadrature
#'
#' @name Scores
NULL



# abstract class structure for scores
setClass("Score", representation(label = "character"))
setClass("ConditionalScore", contains = "Score")
setClass("UnconditionalScore", contains = "Score")

# class for expected scores
setClass("IntegralScore", representation(
        cs                = "ConditionalScore",
        data_distribution = "DataDistribution",
        prior             = "Prior"
    ),
    contains = "UnconditionalScore")



setMethod("show", signature(object = "Score"), function(object) {
    cat(print(object), "\n")
})

setMethod("print", signature('Score'), function(x, ...) {
    if (!is.na(x@label)) x@label else paste0(class(x)[1], "<>")
})

setMethod("print", signature('IntegralScore'), function(x, ...) {
    if (is.na(x@label)) {
        sprintf("E[%s]<%s;%s>", print(x@cs), print(x@data_distribution), print(x@prior))
    } else {
        sprintf("E[%s]<%s;%s>", x@label, print(x@data_distribution), print(x@prior))
    }
})



#' @rdname Scores
#' @export
setGeneric("expected", function(s, data_distribution, prior, ...) standardGeneric("expected"))



#' @rdname Scores
#' @export
setMethod("expected", signature("ConditionalScore"),
          function(s, data_distribution, prior, label = NA_character_, ...) {
              new("IntegralScore", label = label, cs = s,
                  data_distribution = data_distribution, prior = prior)
          })






#' @rdname Scores
#' @export
setGeneric("evaluate", function(s, design, ...) standardGeneric("evaluate"))



#' @rdname Scores
#' @export
setMethod("evaluate", signature("IntegralScore", "TwoStageDesign"),
    function(s, design, optimization = FALSE, subdivisions = 10000L, ...) {
        # use design-specific implementation
        if (optimization) return(.evaluate(s, design, ...))
        # use generic approach
        epsilon <- sqrt(.Machine$double.eps)
        n1      <- design@n1 # n1(design)
        early_futility <- predictive_cdf(s@data_distribution, s@prior, design@c1f, n1) *
            evaluate(s@cs, design, design@c1f - epsilon) # score is constant on early stopping region
        early_efficacy <- (1 - predictive_cdf(s@data_distribution, s@prior, design@c1e, n1)) *
            evaluate(s@cs, design, design@c1e + epsilon)
        continuation <- stats::integrate(
            function(x1) predictive_pdf(s@data_distribution, s@prior, x1, n1, ...) * evaluate(s@cs, design, x1, ...),
            lower        = design@c1f,
            upper        = design@c1e,
            subdivisions = subdivisions
        )$value
        return(early_futility + continuation + early_efficacy)
    })





# not user facing
setGeneric(".evaluate", function(s, design, ...) standardGeneric(".evaluate"))



# not user facing!
setMethod(".evaluate", signature("IntegralScore", "TwoStageDesign"),
    function(s, design, ...) {
        epsilon <- sqrt(.Machine$double.eps)
        n1 <- design@n1; c1f <- design@c1f; c1e <- design@c1e
        early_futility <- predictive_cdf(s@data_distribution, s@prior, c1f, n1) *
            evaluate(s@cs, design, c1f - epsilon, optimization = TRUE, ...) # score is constant on early stopping region
        early_efficacy <- (1 - predictive_cdf(s@data_distribution, s@prior, c1e, n1)) *
            evaluate(s@cs, design, c1e + epsilon, optimization = TRUE, ...)
        continuation <- integrate_rule(
            function(x1) predictive_pdf(s@data_distribution, s@prior, x1, n1, ...) * evaluate(s@cs, design, x1, optimization = TRUE, ...),
            low     = c1f,
            up      = c1e,
            x       = design@x1_norm_pivots,
            weights = design@weights
        )
        return(early_futility + continuation + early_efficacy)
    })



# not user facing!
setMethod(".evaluate", signature("IntegralScore", "OneStageDesign"),
    function(s, design, ...) {
        epsilon <- sqrt(.Machine$double.eps)
        n1 <- design@n1; c1f <- design@c1f; c1e <- design@c1e
        futility <- predictive_cdf(s@data_distribution, s@prior, c1f, n1) *
            evaluate(s@cs, design, c1f - epsilon, optimization = TRUE, ...) # score is constant on early stopping region
        efficacy <- (1 - predictive_cdf(s@data_distribution, s@prior, c1e, n1)) *
            evaluate(s@cs, design, c1e + epsilon, optimization = TRUE, ...)
        return(futility + efficacy)
    })
