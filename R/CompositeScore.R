setClass("CompositeScore", representation(
        expr       = "{",
        scores     = "list",
        non_scores = "list",
        label      = "character"
    ),
    contains = "Score"
)

setClass("CompositeUnconditionalScore",
     contains = c("UnconditionalScore", "CompositeScore")
)

setClass("CompositeConditionalScore",
     contains = c("ConditionalScore", "CompositeScore")
)


setMethod("print", signature("CompositeScore"), function(x, ...) {
    labels <- lapply(x@scores, function(x) utils::capture.output(show((x))))
    str <- as.character(x@expr)
    for (i in 2:length(str)) {
        str[i] <- sub('\\{', "\\{\\{", str[i])
        str[i] <- sub('\\}', "\\}\\}", str[i])
        for (j in 1:length(labels)) {
            str[i] <- sub(names(labels)[j], sprintf("{%s}", names(labels)[j]), str[i])
        }
    }
    glue::glue(paste0(str[-1], collapse = "; "), .envir = c(labels, x@non_scores))
})


#' Score Composition
#'
#' \code{composite} defines new composite scores by point-wise evaluation of
#' scores in any valid numerical expression.
#'
#' @param expr Expression (in curly brackets); must contain at least one score
#'   variable; if multiple scores are used, they must either all be conditional
#'   or unconditional. Currently, no non-score variables are supported
#' @param s object of class \code{CompositeScore}
#' @template design
#' @template label
#' @template dotdotdot
#'
#' @return an object of class \code{CompositeConditionalScore} or
#'   \code{CompositeUnconditionalScore} depending on the class of the scores used
#'   in \code{expr}
#'
#' @seealso \link{Scores}
#'
#' @examples
#' ess   <- ExpectedSampleSize(Normal(), PointMassPrior(.4, 1))
#' power <- Power(Normal(), PointMassPrior(.4, 1))
#'
#' # linear combination:
#' composite({ess - 50*power})
#'
#' # control flow (e.g. for and while loops)
#' composite({
#'   res <- 0
#'   for (i in 1:3) {
#'      res <- res + ess
#'   }
#'   res
#' })
#'
#' # functional composition
#' composite({log(ess)})
#' cp <- ConditionalPower(Normal(), PointMassPrior(.4, 1))
#' composite({3*cp})
#'
#' @export
composite <- function(expr, label = NA_character_) {

    vars <- mget(
        x          = all.vars(substitute(expr)),
        envir      = parent.frame(n = 1),
        inherits   = TRUE, # make sure to look in entire stack
        ifnotfound = list(NA)
    )
    # extract 'Score' variables
    idx_scores <- as.logical(sapply(vars, function(x) is(x, "Score")))
    scores     <- vars[idx_scores]
    non_scores <- vars[!idx_scores]

    if (length(scores) == 0) {
        stop("no scores in expression")
    } else {
        if (any(sapply(scores, function(x) is(x, "ConditionalScore")))) {
            # if we have any conditional score, all must be!
            if (any(sapply(scores, function(x) !is(x, "ConditionalScore")))) {
                stop("either all or none of the scores must be conditional!")
            }
            return(new("CompositeConditionalScore",
                       expr       = substitute(expr),
                       scores     = scores,
                       non_scores = non_scores,
                       label      = label))
        } else {
            # only unconditional scores
            return(new("CompositeUnconditionalScore",
                       expr       = substitute(expr),
                       scores     = scores,
                       non_scores = non_scores,
                       label      = label))
        }
    }

}


#' @rdname composite
#' @export
setMethod("evaluate", signature("CompositeScore", "TwoStageDesign"),
          function(s, design, ...) {
             values <- lapply(s@scores, function(s) evaluate(s, design, ...))
             return(eval(s@expr, c(values, s@non_scores)))
          })
