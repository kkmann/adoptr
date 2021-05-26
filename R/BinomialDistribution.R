#' Binomial data distribution
#'
#' Implements the normal approximation for a test on rates.
#' The reponse rate in the control group,
#' \ifelse{html}{\out{r<sub>C</sub>}}{\eqn{r_C}}, has to be specified by
#' \code{rate_control}.
#' The null hypothesis is:
#' \ifelse{html}{\out{r<sub>E</sub> &le; r<sub>C</sub>}}{\eqn{r_E <= r_C}},
#' where \ifelse{html}{\out{r<sub>E</sub>}}{\eqn{r_E}} denotes the response rate
#' in the invervention group.
#' It is tested against the alternative
#' \ifelse{html}{\out{r<sub>E</sub> > r<sub>C</sub>}}{\eqn{r_E > r_C}}.
#' The test statistic is given as
#' \ifelse{html}{\out{X<sub>1</sub> = &radic;n (r<sub>E</sub> - r<sub>C</sub>) / &radic;(2  r<sub>0</sub> (1-r<sub>0</sub>))}}{\eqn{X_1 = \sqrt{n}(r_E - r_C) / \sqrt{2 r_0 (1- r_0)}}},
#' where \ifelse{html}{\out{r<sub>0</sub>}}{\eqn{r_0}} denotes the mean between
#' \ifelse{html}{\out{r<sub>E</sub>}}{\eqn{r_E}} and
#' \ifelse{html}{\out{r<sub>C</sub>}}{\eqn{r_C}} in the two-armed case,
#' and \ifelse{html}{\out{r<sub>E</sub>}}{\eqn{r_E}} in the one-armed case.#'
#' All priors have to be defined for the rate difference
#' \ifelse{html}{\out{r<sub>E</sub> - r<sub>C</sub>}}{\eqn{r_E - r_C}}.
#'
#' @slot rate_control cf. parameter 'rate_control'
#'
#' @template DataDistributionTemplate
#'
#' @rdname BinomialDataDistribution-class
#' @exportClass Binomial
setClass("Binomial", representation(
    rate_control = "numeric",
    two_armed    = "logical"
),
contains = "DataDistribution")


#' @param rate_control assumed response rate in control group
#' @param two_armed logical indicating if a two-armed trial is regarded
#'
#' @examples
#' datadist <- Binomial(rate_control = 0.2, two_armed = FALSE)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname BinomialDataDistribution-class
#' @export
Binomial <- function(rate_control, two_armed = TRUE) {
    if (any(rate_control >= 1, rate_control <= 0))
        stop("The response rate in the control group must be in (0,1)!")
    new("Binomial", rate_control = rate_control, two_armed = two_armed)
}


#' @examples
#' probability_density_function(Binomial(.2, FALSE), 1, 50, .3)
#'
#' @details If the distribution is \code{\link{Binomial}},
#'   \ifelse{html}{\out{theta}}{\eqn{theta}} denotes the rate difference between
#'   intervention and control group.
#'   Then, the mean is assumed to be
#'   \ifelse{html}{\out{&radic; n  theta}}{\eqn{\sqrt{n} theta}}.
#'
#' @rdname probability_density_function
#' @export
setMethod("probability_density_function", signature("Binomial", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              rate_intervention <- theta + dist@rate_control
              if (any(rate_intervention >= 1, rate_intervention <= 0)) # nocov start
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              # nocov end

              sigma_A <- sqrt(rate_intervention * (1 - rate_intervention) +
                                  ifelse(dist@two_armed, dist@rate_control * (1 - dist@rate_control), 0))
              p_0     <- (rate_intervention + ifelse(dist@two_armed, dist@rate_control, rate_intervention)) / 2
              sigma_0 <- sqrt(2 * p_0 * (1 - p_0))

              return(stats::dnorm(x, mean = sqrt(n) * theta / sigma_0, sd = sigma_A / sigma_0))
          })


#' @examples
#' cumulative_distribution_function(Binomial(.1, TRUE), 1, 50, .3)
#'
#' @details If the distribution is \code{\link{Binomial}},
#'   \ifelse{html}{\out{theta}}{\eqn{theta}} denotes the rate difference between
#'   intervention and control group.
#'   Then, the mean is assumed to be
#'   \ifelse{html}{\out{&radic; n  theta}}{\eqn{\sqrt{n} theta}}.
#'
#' @rdname cumulative_distribution_function
#' @export
setMethod("cumulative_distribution_function", signature("Binomial", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
              rate_intervention <- theta + dist@rate_control
              if (any(rate_intervention >= 1, rate_intervention <= 0)) # nocov start
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              # nocov end
              sigma_A <- sqrt(rate_intervention * (1 - rate_intervention) +
                                  ifelse(dist@two_armed, dist@rate_control * (1 - dist@rate_control), 0))
              p_0     <- (rate_intervention + ifelse(dist@two_armed, dist@rate_control, rate_intervention)) / 2
              sigma_0 <- sqrt(2 * p_0 * (1 - p_0))

              return(stats::pnorm(x, mean = sqrt(n) * theta / sigma_0, sd = sigma_A / sigma_0))
        })



#' @param probs vector of probabilities
#' @rdname BinomialDataDistribution-class
#' @export
setMethod("quantile", signature("Binomial"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
              rate_intervention <- theta + x@rate_control
              if (any(rate_intervention >= 1, rate_intervention <= 0)) # nocov start
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              # nocov end
              sigma_A <- sqrt(rate_intervention * (1 - rate_intervention) +
                                  ifelse(x@two_armed, x@rate_control * (1 - x@rate_control), 0))
              p_0     <- (rate_intervention + ifelse(x@two_armed, x@rate_control, rate_intervention)) / 2
              sigma_0 <- sqrt(2 * p_0 * (1 - p_0))

              return(stats::qnorm(probs, mean = sqrt(n) * theta / sigma_0, sd = sigma_A / sigma_0))
          })



#' @details Note that \code{simulate} for class \code{Binomial} simulates the
#'    normal approximation of the test statistic.
#'
#' @rdname BinomialDataDistribution-class
#'
#' @param object object of class \code{Binomial}
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("Binomial", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
              rate_intervention <- theta + object@rate_control
              if (any(rate_intervention >= 1, rate_intervention <= 0)) # nocov start
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              # nocov end
              sigma_A <- sqrt(rate_intervention * (1 - rate_intervention) +
                                  ifelse(object@two_armed, object@rate_control * (1 - object@rate_control), 0))
              p_0     <- (rate_intervention + ifelse(object@two_armed, object@rate_control, rate_intervention)) / 2
              sigma_0 <- sqrt(2 * p_0 * (1 - p_0))

              if (!is.null(seed)) set.seed(seed)

              return(stats::rnorm(nsim, mean = sqrt(n) * theta / sigma_0, sd = sigma_A / sigma_0))
})




setMethod("print", signature('Binomial'), function(x, ...) {
    glue::glue(
        "{class(x)[1]}<{if (x@two_armed) 'two-armed' else 'single-armed'}>",
        'response rate in control group: {x@rate_control}',
        .sep = ", "
    )
})
