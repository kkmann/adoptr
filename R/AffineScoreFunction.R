setClass("AffineScore", representation(
    scores    = "list",
    coefs     = "numeric",
    intercept = "numeric"
    ))


AffineScore <- function(scores, coefs, intercept) {
    if (length(scores) != length(coefs))
        stop("scores and coefs must have same size")
    if (length(intercept) != 1)
        stop("intercept must have length 0")
    if (any(!is.finite(c(coefs, intercept))))
        stop("scores and intercept must be finite")
    new("AffineScore", scores = scores, coefs = coefs, intercept = intercept)
}


setMethod("evaluate", signature("AffineScore", "Design"),
          function(s, design, ...) sum(s@coefs * sapply(s@scores, function(s, ...) evaluate(s, design, ...), ...) + s@intercept) )




setClass("AffineUnconditionalScore", contains = c("AffineScore", "UnconditionalScore"))

AffineUnconditionalScore <- function(scores, coefs, intercept = 0) {
    if (!all(sapply(scores, function(s) is(s, "UnconditionalScore"))))
        stop("all scores must be unconditional scores")
    res <- AffineScore(scores, coefs, intercept)
    class(res) <- "AffineUnconditionalScore"
    return(res)
}

setMethod("+", signature("AffineUnconditionalScore", "UnconditionalScore"),
          function(e1, e2) AffineUnconditionalScore(c(e1@scores, list(e2)), c(e1@coefs, 1.0), e1@intercept) )
setMethod("+", signature("UnconditionalScore", "AffineUnconditionalScore"),
          function(e1, e2) e2 + e1 )

setMethod("+", signature("AffineUnconditionalScore", "numeric"),
          function(e1, e2) AffineUnconditionalScore(e1@scores, e1@coefs, e1@intercept + e2) )
setMethod("+", signature("numeric", "AffineUnconditionalScore"),
          function(e1, e2) e2 + e1 )

## TODO: check for duplicate scores and combine!
setMethod("+", signature("AffineUnconditionalScore", "AffineUnconditionalScore"),
          function(e1, e2) AffineUnconditionalScore(c(e1@scores, e2@scores), c(e1@coefs, e2@coefs), e1@intercept + e2@intercept) )

setMethod("*", signature("AffineUnconditionalScore", "numeric"),
          function(e1, e2) AffineUnconditionalScore(e1@scores, e2 * e1@coefs, e2 * e1@intercept) )
setMethod("*", signature("numeric", "AffineUnconditionalScore"),
          function(e1, e2) e2 * e1 )



setClass("AffineConditionalScore", contains = c("AffineScore", "AbstractConditionalScore"))

AffineConditionalScore <- function(scores, coefs, intercept = 0) {
    if (!all(sapply(scores, function(s) is(s, "ConditionalScore"))))
        stop("all scores must be conditional scores")
    res <- AffineScore(scores, coefs, intercept)
    class(res) <- "AffineConditionalScore"
    return(res)
}

setMethod("+", signature("AffineConditionalScore", "ConditionalScore"),
          function(e1, e2) AffineConditionalScore(c(e1@scores, list(e2)), c(e1@coefs, 1.0), e1@intercept) )
setMethod("+", signature("ConditionalScore", "AffineConditionalScore"),
          function(e1, e2) e2 + e1 )

setMethod("+", signature("AffineConditionalScore", "numeric"),
          function(e1, e2) AffineConditionalScore(e1@scores, e1@coefs, e1@intercept + e2) )
setMethod("+", signature("numeric", "AffineConditionalScore"),
          function(e1, e2) e2 + e1 )

## TODO: check for duplicate scores and combine!
setMethod("+", signature("AffineConditionalScore", "AffineConditionalScore"),
          function(e1, e2) AffineConditionalScore(c(e1@scores, e2@scores), c(e1@coefs, e2@coefs), e1@intercept + e2@intercept) )

setMethod("*", signature("AffineConditionalScore", "numeric"),
          function(e1, e2) AffineConditionalScore(e1@scores, e2 * e1@coefs, e2 * e1@intercept) )
setMethod("*", signature("numeric", "AffineConditionalScore"),
          function(e1, e2) e2 * e1 )
