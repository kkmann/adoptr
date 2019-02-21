#' Conditional power of a design given stage-one outcome
#'
#' This score evaluates to \eqn{P[X_2 > c2(design, X_1)|X_1 = x_1]}.
#' Note that the distribution of \eqn{X_2} is the posterior predictive after
#' observing \eqn{X_1=x_1}.
#'
#' @details See documentation in \link{ConditionalScore-class}
#'     for further details and inherited methods.
#'
#' @template ConditionalScoreTemplate
#' @param dist data distribution
#' @param prior prior distribution
#'
#' @exportClass ConditionalPower
setClass("ConditionalPower", contains = "ConditionalScore")

#' @describeIn ConditionalPower-class constructor
#' @export
ConditionalPower <- function(dist, prior) new("ConditionalPower", distribution = dist, prior = prior)

#' @rdname ConditionalPower-class
setMethod("evaluate", signature("ConditionalPower", "TwoStageDesign"),
          function(s, design, x1, optimization = FALSE, ...) {
              sapply(x1,
                  function(x1) expectation(
                      posterior(s@distribution, s@prior, x1, design@n1, ...),
                      function(theta)
                          1 - cumulative_distribution_function(s@distribution, c2(design, x1), n2(design, x1), theta)
                  )
              )
          })
