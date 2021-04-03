#' Normal data distribution
#'
#' Implements a normal data distribution for z-values given an observed z-value
#' and stage size.
#' Standard deviation is 1 and mean \ifelse{html}{\out{&theta; &radic;n}}{\eqn{\theta\sqrt n}} where
#' \ifelse{html}{\out{&theta;}}{\eqn{\theta}} is the standardized effect size.
#' The option \code{two_armed} can be set to decide whether a one-arm or a
#' two-arm design should be computed.
#'
#' See \code{\link{DataDistribution-class}} for more details.
#'
#' @template DataDistributionTemplate
#'
#' @rdname NormalDataDistribution-class
#' @exportClass Normal
setClass("Normal", representation(
    two_armed = "logical"
),
contains = "DataDistribution")


#' @param two_armed logical indicating if a two-armed trial is regarded
#'
#' @examples
#' datadist <- Normal(two_armed = TRUE)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname NormalDataDistribution-class
#' @export
Normal <- function(two_armed = TRUE) new("Normal", two_armed = two_armed)


#' @examples
#' probability_density_function(Normal(), 1, 50, .3)
#'
#' @details If the distribution is \code{\link{Normal}}, then
#'   the mean is assumed to be
#'   \ifelse{html}{\out{&radic; n  theta}}{\eqn{\sqrt{n} theta}}.
#'
#' @rdname probability_density_function
#' @export
setMethod("probability_density_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              if (dist@two_armed) {
                  theta <- theta / sqrt(2)
              }
              stats::dnorm(x, mean = sqrt(n) * theta, sd = 1)
          })


#' @examples
#' cumulative_distribution_function(Normal(), 1, 50, .3)
#'
#' @details If the distribution is \code{\link{Normal}}, then
#'   the mean is assumed to be
#'   \ifelse{html}{\out{&radic; n  theta}}{\eqn{\sqrt{n} theta}}.
#'
#' @rdname cumulative_distribution_function
#' @export
setMethod("cumulative_distribution_function", signature("Normal", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              if (dist@two_armed) {
                  theta <- theta / sqrt(2)
              }
              stats::pnorm(x, mean = sqrt(n) * theta, sd = 1)
          })



#' @param probs vector of probabilities
#' @rdname NormalDataDistribution-class
#' @export
setMethod("quantile", signature("Normal"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
              if (x@two_armed) {
                  theta <- theta / sqrt(2)
              }
              stats::qnorm(probs, mean = sqrt(n) * theta, sd = 1)
          })



#' @rdname NormalDataDistribution-class
#'
#' @param object object of class \code{Normal}
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("Normal", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
              if (object@two_armed)
                  theta <- theta / sqrt(2)

              if (!is.null(seed))
                  set.seed(seed)

              stats::rnorm(nsim, mean = sqrt(n) * theta, sd = 1)
          })



setMethod("print", signature('Normal'), function(x, ...) {
    glue::glue(
        "{class(x)[1]}<{if (x@two_armed) 'two-armed' else 'single-armed'}>"
    )
})
