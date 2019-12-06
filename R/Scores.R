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



# internal method to evaluate condition scores effectively on grid
setGeneric(".evaluate", function(s, design, ...) standardGeneric(".evaluate"))

setMethod(".evaluate", signature("ConditionalScore", "TwoStageDesign"),
          function(s, design, ...) {
              epsilon <- sqrt(.Machine$double.eps)
              pivots  <- scaled_integration_pivots(design)
              return(list(
                  early_futility = evaluate(s, design, design@c1f - epsilon, optimization = TRUE, ...),
                  early_efficacy = evaluate(s, design, design@c1e + epsilon, optimization = TRUE, ...),
                  pivots         = evaluate(s, design, pivots, optimization = TRUE, ...)
              ))
            })




#' @rdname Scores
#' @export
setGeneric("evaluate", function(s, design, ...) standardGeneric("evaluate"))



#' @rdname Scores
#' @export
setMethod("evaluate", signature("IntegralScore", "TwoStageDesign"),
    function(s, design, optimization = FALSE, subdivisions = 10000L, ...) {
        epsilon <- sqrt(.Machine$double.eps)
        c1f     <- design@c1f; c1e <- design@c1e
        n1      <- if (optimization) design@n1 else n1(design, round = TRUE)
        pr_ef   <- predictive_cdf(s@data_distribution, s@prior, c1f, n1)
        pr_ee   <- 1 - predictive_cdf(s@data_distribution, s@prior, c1e, n1)
        if (optimization == TRUE) {
            cs             <- .evaluate(s@cs, design)
            early_stopping <- pr_ef * cs$early_futility + pr_ee * cs$early_efficacy
            if (is(design, 'OneStageDesign')) return(early_stopping)
            pivots         <- scaled_integration_pivots(design)
            pdf            <- predictive_pdf(s@data_distribution, s@prior, pivots, n1)
            continuation   <- gauss_quad(cs$pivots * pdf, c1f, c1e, design@weights)
        } else {
            early_stopping <- pr_ef * evaluate(s@cs, design, c1f - epsilon, optimization = FALSE, ...) +
                pr_ee * evaluate(s@cs, design, c1e + epsilon, optimization = FALSE, ...)
            if (is(design, 'OneStageDesign')) return(early_stopping)
            continuation <- stats::integrate(
                function(x1) predictive_pdf(s@data_distribution, s@prior, x1, n1, ...) *
                    evaluate(s@cs, design, x1, ...),
                lower = c1f, upper = c1e, subdivisions = subdivisions
            )$value
        }
        return(early_stopping + continuation)
    })
