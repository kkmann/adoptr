#' Log-rank test
#'
#' Implements the normal approximation of the log-rank test statistic.
#'
#' @slot event_rate cf. parameter 'event_rate'
#' @template DataDistributionTemplate
#'
#' @rdname SurvivalDataDistribution-class
#' @exportClass Survival
setClass("Survival", representation (
    event_rate = "numeric",
    two_armed = "logical"),contains = "DataDistribution")



#' @param event_rate probability that a subject will eventually have an event
#' @param two_armed logical indicating if a two-armed trial is regarded
#'
#' @examples
#' datadist <- Survival(event_rate=0.6, two_armed=TRUE)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname SurvivalDataDistribution-class
#' @export
Survival <- function(event_rate,two_armed=TRUE){
    if(any(event_rate>=1,event_rate<=0)){
    stop("The assumed event rate must be in (0,1)!")}
    new("Survival",event_rate = event_rate,two_armed = two_armed)
}

#' @examples
#' probability_density_function(Survival(0.6,TRUE),0.75,50,0.9)
#'
#' @rdname probability_density_function
#' @export
setMethod("probability_density_function",
          signature("Survival","numeric","numeric","numeric"),
          function(dist,x,n,theta, ...){
              if(dist@two_armed) n <- n/2
              return(stats::dnorm(x,mean = sqrt(n)*log(theta),sd = 1) )
})

#' @examples
#' cumulative_distribution_function(Survival(0.6,TRUE),0.75,50,0.9)
#'
#' @rdname cumulative_distribution_function
#' @export
setMethod("cumulative_distribution_function",
          signature("Survival","numeric","numeric","numeric"),
          function(dist, x, n, theta, ...){
              if(dist@two_armed) n <- n/2
              return(stats::pnorm(x,mean=sqrt(n)*log(theta),sd=1))
})

#' @param probs vector of probabilities
#' @rdname SurvivalDataDistribution-class
#' @export
setMethod("quantile", signature("Survival"),
          function(x, probs, n, theta, ...) {
              if(x@two_armed) n <- n/2
              return(stats::qnorm(probs, mean=sqrt(n)*log(theta),sd=1))
})

#' @rdname SurvivalDataDistribution-class
#'
#' @param object object of class \code{Survival}
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("Survival", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
              if(object@two_armed) n <- n/2
              if (!is.null(seed)) set.seed(seed)
              return(stats::rnorm(nsim, mean=sqrt(n)*log(theta), sd = 1))
          })

setMethod("print", signature('Survival'), function(x, ...) {
    glue::glue(
        "{class(x)[1]}<{if (x@two_armed) 'two-armed' else 'single-armed'}>",
        "event rate: {x@event_rate}",
        .sep=", "
    )
})
