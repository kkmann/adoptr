#' @export
setClass("CompositeScore", list(
    expr = "{",
    vars = "list"
),
contains = "Score"
)

#' @export
CompositeScore <- function(expr) {

    vars <- mget(
        x          = all.vars(substitute(expr)),
        envir      = parent.frame(n = 1),
        inherits   = FALSE,
        ifnotfound = list(NA)
    )

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
