setClass("CompositeScore", representation(
        expr   = "{",
        scores = "list"
    ),
    contains = "Score"
)

setClass("CompositeUnconditionalScore",
     contains = c("UnconditionalScore", "CompositeScore")
)

setClass("CompositeConditionalScore",
     contains = c("ConditionalScore", "CompositeScore")
)

#' @export
compose <- function(expr) {

    vars <- mget(
        x          = all.vars(substitute(expr)),
        envir      = parent.frame(n = 1),
        inherits   = TRUE, # make sure to look in entire stack
        ifnotfound = list(NA)
    )
    # extract 'Score' variables
    scores <- vars[as.logical(sapply(vars, function(x) is(x, "Score")))]

    if (length(scores) == 0) {
        stop("no scores in expression")
    } else {
        if (any(sapply(scores, function(x) is(x, "ConditionalScore")))) {
            # if we have any conditional score, all must be!
            if (any(sapply(scores, function(x) !is(x, "ConditionalScore")))) {
                stop("either all or none of the scores must be conditional!")
            }
            return(new("CompositeConditionalScore",
                       expr         = substitute(expr),
                       scores       = scores
            ))
        } else {
            # only unconditional scores
            return(new("CompositeUnconditionalScore", expr = substitute(expr), scores = scores))
        }
    }

}



#' @export
setMethod("evaluate", signature("CompositeScore", "TwoStageDesign"),
          function(s, design, ...) {
             values <- lapply(s@scores, function(s) evaluate(s, design, ...))
             return(eval(s@expr, values))

          })
