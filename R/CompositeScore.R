#' @export
setClass("CompositeScore", list(
    expr = "{",
    vars = "list"
),
contains = "Score"
)

#' @export
compose <- function(expr) {

    vars <- mget(
        x          = all.vars(substitute(expr)),
        envir      = parent.frame(n = 1),
        inherits   = FALSE,
        ifnotfound = list(NA)
    )

    # check if ConditionalScores are in vars (do we still need AbstractConditionalScores?)
    # if yes:
    #   all scores must be ConditionalScores with same prior / data distribution
    #   return CompositeConditionalScore (subclass of ConditionalScore
    #   -> that's why we need unique prior/data distribution)
    # if no:
    #   only unconditional scores, ok return CompositeScore (subsclass of Score?)

    res <- new("CompositeScore",
               expr  = substitute(expr),
               vars  = vars
    )

    return(res)

}


#' @export
setMethod("evaluate", signature("CompositeScore", "TwoStageDesign"),
          function(s, design, ...) {

              values = list()
              if (length(s@vars) > 0) {
                  values   <- lapply(
                      s@vars[sapply(s@vars, function(x) is(x, "Score"))],
                      function(s) evaluate(s, design, ...)
                  )
              }

              return(eval(s@expr, values))

          }
)
