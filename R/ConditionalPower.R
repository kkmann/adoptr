#' Conditional power of a design given stage-one outcome
#'
#' This score evaluates \ifelse{html}{\out{P[X<sub>2</sub> > c2(design, X<sub>1</sub>) | X<sub>1</sub> = x<sub>1</sub>]}}{\eqn{\boldsymbol{P}[X_2 > c_2(design, X_1)|X_1 = x_1]}}.
#' Note that the distribution of \ifelse{html}{\out{X<sub>2</sub>}}{\eqn{X_2}} is the posterior predictive after
#' observing \ifelse{html}{\out{X<sub>1</sub> = x<sub>1</sub>}}{\eqn{X_1 = x_1}}.
#'
#'
#' @template dist
#' @template prior
#'
#' @seealso The method \code{\link{evaluate}} provides evaluation of the
#'    \code{ConditionalPower}.
#'
#' @aliases ConditionalPower
#' @exportClass ConditionalPower
setClass("ConditionalPower", contains = "ConditionalScore")

#' @examples
#' cp <- ConditionalPower(dist = Normal(), prior = PointMassPrior(.4, 1))
#'
#' @rdname ConditionalPower-class
#' @export
ConditionalPower <- function(dist, prior) new("ConditionalPower", distribution = dist, prior = prior)

#' @examples
#' # evaluate conditional power
#' evaluate(
#'    ConditionalPower(Normal(), PointMassPrior(.3, 1)),
#'    TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L),
#'    x1 = 1
#' )
#'
#' @template x1
#'
#' @rdname evaluate
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
