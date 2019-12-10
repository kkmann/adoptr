#' Data distributions
#'
#' \code{DataDistribution} is an abstract class used to represent the distribution
#' of a sufficient statistic \code{x} given a sample size \code{n} and a
#' single parameter value \code{theta}.
#'
#' This abstraction layer allows the representation of t-distributions
#' (unknown variance), normal distribution (known variance), and normal
#' approximation of a binary endpoint.
#' Currently, only the normal case is implemented with \code{\link{Normal-class}}.
#'  TODO add binomial
#'
#' The logical option \code{two_armed} allows to decide whether a one-arm or
#' a two-arm (the default) design should be computed. In the case of a two-arm
#' design all sample sizes are per group.
#'
#' @slot two_armed Logical that indicates if a two-arm design is assumed.
#'
#' @examples
#' datadist <- Normal(two_armed = FALSE)
#'
#' @template DataDistributionTemplate
#'
#' @aliases DataDistribution
#' @exportClass DataDistribution
setClass("DataDistribution", representation(
    two_armed = "logical")
)


#' Probability density function
#'
#' \code{probability_density_function} evaluates the probability density
#' function of a specific distribution \code{dist} at a point \code{x}.
#'
#' @template dist
#' @template DataDistributionTemplate
#'
#' @export
setGeneric("probability_density_function", function(dist, x, n, theta, ...) standardGeneric("probability_density_function"))


#' Cumulative distribution function
#'
#' \code{cumulative_distribution_function} evaluates the cumulative distribution
#' function of a specific distribution \code{dist} at a point \code{x}.
#'
#' @template dist
#' @template DataDistributionTemplate
#'
#' @export
setGeneric("cumulative_distribution_function", function(dist, x, n, theta, ...) standardGeneric("cumulative_distribution_function"))





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

setMethod("show", signature(object = "DataDistribution"), function(object) {
    cat(print(object), "\n")
})





#' Binomial data distribution
#'
#' Implements the normal approximation for a test on rates.
#' The reponse rate in the control group,
#' \ifelse{html}{\out{r<sub>C</sub>}}{\eqn{r_C}}, has to be specified by
#' \code{rate_control}.
#' The null hypothesis is:
#' \ifelse{html}{\out{r<sub>E</sub> &le; r<sub>E</sub>}}{\eqn{r_E <= r_C}},
#' where \ifelse{html}{\out{r<sub>E</sub>}}{\eqn{r_E}} denotes the response rate
#' in the invervention group.
#' It is tested against the alternative
#' \ifelse{html}{\out{r<sub>E</sub> > r<sub>E</sub>}}{\eqn{r_E > r_C}}
#'
#' All priors have to be defined for the rate difference
#' \ifelse{html}{\out{r<sub>E</sub> - r<sub>C</sub>}}{\eqn{r_E - r_C}}.
#'
#' \code{Binomial} contains the standard deviation under the null hypothesis
#' to save runtime.
#'
#' @slot rate_control assumed response rate in control group
#' @slot sigma_0 standard deviation of rate difference under the null hypothesis
#'
#' @template DataDistributionTemplate
#'
#' @rdname BinomialDataDistribution-class
#' @exportClass Binomial
setClass("Binomial", representation(
    rate_control = "numeric",
    sigma_0      = "numeric",
    two_armed    = "logical"
),
contains = "DataDistribution")


#' @param rate_control cf. slot 'rate_control'
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
    if(rate_control >= 1 || rate_control <= 0)
        stop("The response rate in the control group must be in (0,1)!")
    sigma_0 = 2 * sqrt(rate_control * (1 - rate_control))
    sigma_0 = ifelse(two_armed, sigma_0, sigma_0 / 2)
    new("Binomial", rate_control = rate_control, two_armed = two_armed, sigma_0 = sigma_0)
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
              if(rate_intervention >= 1 || rate_intervention <= 0)
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              sigma_A <- ifelse(dist@two_armed,
                                sqrt(2 * (dist@rate_control * (1 - dist@rate_control) + rate_intervention * (1 - rate_intervention))),
                                sqrt(rate_intervention * (1 - rate_intervention)))
              return(stats::dnorm(x, mean = sqrt(n) * theta / dist@sigma_0, sd = sigma_A / dist@sigma_0))
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
              if(rate_intervention >= 1 || rate_intervention <= 0)
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              sigma_A <- ifelse(dist@two_armed,
                                sqrt(2 * (dist@rate_control * (1 - dist@rate_control) + rate_intervention * (1 - rate_intervention))),
                                sqrt(rate_intervention * (1 - rate_intervention)))
              return(stats::pnorm(x, mean = sqrt(n) * theta / dist@sigma_0, sd = sigma_A / dist@sigma_0))
          })



#' @param probs vector of probabilities
#' @rdname BinomialDataDistribution-class
#' @export
setMethod("quantile", signature("Binomial"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
              rate_intervention <- theta + dist@rate_control
              if(rate_intervention >= 1 || rate_intervention <= 0)
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              sigma_A <- ifelse(dist@two_armed,
                                sqrt(2 * (dist@rate_control * (1 - dist@rate_control) + rate_intervention * (1 - rate_intervention))),
                                sqrt(rate_intervention * (1 - rate_intervention)))
              return(stats::qnorm(probs, mean = sqrt(n) * theta / dist@sigma_0, sd = sigma_A / dist@sigma_0))
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
              if(rate_intervention >= 1 || rate_intervention <= 0)
                  stop("The response rate in the intervention group must be in (0,1)! Probably the combination of prior and control rate is ill-defined.")
              sigma_A <- ifelse(object@two_armed,
                                sqrt(2 * (object@rate_control * (1 - object@rate_control) + rate_intervention * (1 - rate_intervention))),
                                sqrt(rate_intervention * (1 - rate_intervention)))

              if (!is.null(seed))
                  set.seed(seed)

              stats::rnorm(nsim, mean = sqrt(n) * theta / object@sigma_0, sd = sigma_A / object@sigma_0)
})




setMethod("print", signature('Binomial'), function(x, ...) {
    glue::glue(
        "{class(x)[1]}<{if (x@two_armed) 'two-armed' else 'single-armed'}>",
        'response rate in control group: {x@rate_control}',
        .sep = ", "
    )
})
