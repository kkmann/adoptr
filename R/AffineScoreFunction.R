#' Affine functions of scores
#'
#' These methods allow to add and multiplicate two or several
#' conditional resp. unconditional scores.
#' They can also be weighted by numeric coefficients and shifted
#' by a numeric intercept.
#' This yields to new scores as sum of scores and is therefore a useful
#' tool to combine scores.
#'
#' @slot scores the list of scores
#' @slot coefs numeric vector of the same length as \code{scores}, holding the
#'     coefficients
#' @slot intercept the intercept for the affine function
#'
#' \code{AffineScore} allows to add scores of arbitrary class and shift
#' them by the \code{intercept}.
#' \code{AffineUnconditionalScore} allows the same for scores of class
#' \code{UnconditionalScore} and \code{AffineConditionalScore} for
#' scores of class \code{ConditionalScore}.
#' By the methods \code{+} and \code{*} scores can be added and multiplicated,
#' respectively.
#' Note that it is not possible to use these methods for a mixture of
#' \code{ConditionalScore} and \code{UnconditionalScore} as these
#' require different evaluation techniques.
#' However, for both score classes multiplication and addition with
#' numerics is provided.
#'
#' @exportClass AffineScore
setClass("AffineScore", representation(
    scores    = "list",
    coefs     = "numeric",
    intercept = "numeric"
    ))


#' @param scores cf. corresponding slot
#' @param coefs cf. corresponding slot
#' @param intercept cf. corresponding slot
#' @template dotdotdotTemplate
#'
#' @rdname AffineScore-class
#' @export
AffineScore <- function(scores, coefs, intercept) {
    if (length(scores) != length(coefs))
        stop("scores and coefs must have same size")
    if (length(intercept) != 1)
        stop("intercept must have length 1")
    if (any(!is.finite(c(coefs, intercept))))
        stop("scores and intercept must be finite")
    new("AffineScore", scores = c(scores), coefs = coefs, intercept = intercept)
}



#' @param s score
#' @param design design
#'
#' @rdname AffineScore-class
#' @export
setMethod("evaluate", signature("AffineScore", "TwoStageDesign"),
          function(s, design, ...) {
              res <- 0
              for (i in 1:length(s@scores)) {
                  res <- res + s@coefs[[i]] * evaluate(s@scores[[i]], design, ...) # score might be evaluated at more than one point
              }
              return(res + s@intercept)
          })



#' @rdname AffineScore-class
#'
#' @param object object of class \code{AffineScore}
#' @export
setMethod("show", signature(object = "AffineScore"),
          function(object) cat(class(object)[1]))



setClass("AffineUnconditionalScore", contains = c("AffineScore", "UnconditionalScore"))

AffineUnconditionalScore <- function(scores, coefs, intercept = 0) {
    if (!all(sapply(c(scores), function(s) is(s, "UnconditionalScore"))))
        stop("all scores must be unconditional scores")
    res <- AffineScore(c(scores), coefs, intercept)
    class(res) <- "AffineUnconditionalScore"
    return(res)
}


#'@rdname score-arithmetic
setMethod("+", signature("AffineUnconditionalScore", "UnconditionalScore"),
          function(e1, e2) AffineUnconditionalScore(c(e1@scores, list(e2)), c(e1@coefs, 1.0), e1@intercept) )
#'@rdname score-arithmetic
setMethod("+", signature("UnconditionalScore", "AffineUnconditionalScore"),
          function(e1, e2) e2 + e1 )
#'@rdname score-arithmetic
setMethod("+", signature("AffineUnconditionalScore", "numeric"),
          function(e1, e2) AffineUnconditionalScore(e1@scores, e1@coefs, e1@intercept + e2) )
#'@rdname score-arithmetic
setMethod("+", signature("numeric", "AffineUnconditionalScore"),
          function(e1, e2) e2 + e1 )
## TODO: check for duplicate scores and combine!
#'@rdname score-arithmetic
setMethod("+", signature("AffineUnconditionalScore", "AffineUnconditionalScore"),
          function(e1, e2) AffineUnconditionalScore(c(e1@scores, e2@scores), c(e1@coefs, e2@coefs), e1@intercept + e2@intercept) )


#'@rdname score-arithmetic
setMethod("*", signature("AffineUnconditionalScore", "numeric"),
          function(e1, e2) AffineUnconditionalScore(e1@scores, e2 * e1@coefs, e2 * e1@intercept) )
#'@rdname score-arithmetic
setMethod("*", signature("numeric", "AffineUnconditionalScore"),
          function(e1, e2) e2 * e1 )



setClass("AffineConditionalScore", contains = c("AffineScore", "AbstractConditionalScore"))

AffineConditionalScore <- function(scores, coefs, intercept = 0) {
    if (!all(sapply(c(scores), function(s) is(s, "AbstractConditionalScore"))))
        stop("all scores must be conditional scores")
    res <- AffineScore(c(scores), coefs, intercept)
    class(res) <- "AffineConditionalScore"
    return(res)
}



#'@rdname score-arithmetic
setMethod("+", signature("AffineConditionalScore", "ConditionalScore"),
          function(e1, e2) AffineConditionalScore(c(e1@scores, list(e2)), c(e1@coefs, 1.0), e1@intercept) )
#'@rdname score-arithmetic
setMethod("+", signature("ConditionalScore", "AffineConditionalScore"),
          function(e1, e2) e2 + e1 )
#'@rdname score-arithmetic
setMethod("+", signature("AffineConditionalScore", "numeric"),
          function(e1, e2) AffineConditionalScore(e1@scores, e1@coefs, e1@intercept + e2) )
#'@rdname score-arithmetic
setMethod("+", signature("numeric", "AffineConditionalScore"),
          function(e1, e2) e2 + e1 )
## TODO: check for duplicate scores and combine!
#'@rdname score-arithmetic
setMethod("+", signature("AffineConditionalScore", "AffineConditionalScore"),
          function(e1, e2) AffineConditionalScore(c(e1@scores, e2@scores), c(e1@coefs, e2@coefs), e1@intercept + e2@intercept) )


#'@rdname score-arithmetic
setMethod("*", signature("AffineConditionalScore", "numeric"),
          function(e1, e2) AffineConditionalScore(e1@scores, e2 * e1@coefs, e2 * e1@intercept) )
#'@rdname score-arithmetic
setMethod("*", signature("numeric", "AffineConditionalScore"),
          function(e1, e2) e2 * e1 )
