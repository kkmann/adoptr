#' (Conditional) Power of a Design
#'
#' This score evaluates \ifelse{html}{\out{P[X<sub>2</sub> > c2(design, X<sub>1</sub>) | X<sub>1</sub> = x<sub>1</sub>]}}{\eqn{\boldsymbol{P}[X_2 > c_2(design, X_1)|X_1 = x_1]}}.
#' Note that the distribution of \ifelse{html}{\out{X<sub>2</sub>}}{\eqn{X_2}} is the posterior predictive after
#' observing \ifelse{html}{\out{X<sub>1</sub> = x<sub>1</sub>}}{\eqn{X_1 = x_1}}.
#'
#' @template dist
#' @template prior
#' @template design
#' @template s
#' @template x1
#' @template optimization
#' @template dotdotdot
#'
#' @seealso \code{\link{Scores}}
#'
#' @examples
#' prior <- PointMassPrior(.4, 1)
#' cp <- ConditionalPower(Normal(), prior)
#' evaluate(
#'    cp,
#'    TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L),
#'    x1 = 1
#' )
#' # these two are equivalent:
#' expected(cp, Normal(), prior)
#' Power(Normal(), prior)
#'
#' @aliases ConditionalPower
#' @exportClass ConditionalPower
setClass("ConditionalPower", representation(
        distribution = "DataDistribution",
        prior        = "Prior"
    ),
    contains = "ConditionalScore")


#' @template label
#'
#' @rdname ConditionalPower-class
#' @export
ConditionalPower <- function(dist, prior, label = "Pr[x2>=c2(x1)|x1]") {
    new("ConditionalPower", distribution = dist, prior = prior, label = label)
}



setMethod("print", signature('ConditionalPower'), function(x, ...) {
    name <- if (!is.na(x@label)) x@label else class(x)[1]
    return(sprintf("%s<%s;%s>", name, print(x@distribution), print(x@prior)))
})



#' @rdname ConditionalPower-class
#' @export
Power <- function(dist, prior, label = "Pr[x2>=c2(x1)]") {
    expected(ConditionalPower(dist, prior), dist, prior, label = label)
}



#' @rdname ConditionalPower-class
#' @export
setMethod("evaluate", signature("ConditionalPower", "TwoStageDesign"),
          function(s, design, x1, optimization = FALSE, ...) {
              sapply(x1,
                  function(x1) expectation(
                      posterior(s@distribution, s@prior, x1, n1(design, round = !optimization), ...),
                      function(theta)
                          1 - cumulative_distribution_function(s@distribution, c2(design, x1), n2(design, x1, round = !optimization), theta)
                  )
              )
          })
