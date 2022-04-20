#' Log-rank test
#'
#' Implements the normal approximation of the log-rank test statistic.
#' @template DataDistributionTemplate
#'
#' @rdname LogRank-class
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
#' @rdname LogRank-class
#' @export
Survival <- function(event_rate,two_armed=TRUE){
    if(any(event_rate>=1,event_rate<=0))
    stop("The assumed event rate must be in (0,1)!")
    new("Survival",event_rate = event_rate,two_armed = two_armed)
}
#' @rdname LogRank-class
#' @export
setMethod("probability_density_function",
            signature("Survival","numeric","numeric","numeric"),
            function(dist,x,n,theta,...){
                if(dist@two_armed) n <- n/2
                return(stats::dnorm(x,mean = sqrt(n*dist@event_rate)*log(theta),sd = 1) )
})


#' @rdname LogRank-class
#' @export
setMethod("cumulative_distribution_function",
          signature("Survival","numeric","numeric","numeric"),
          function(dist, x, n, theta, ...){
              if(dist@two_armed) n <- n/2
              return(stats::pnorm(x,mean=sqrt(n*dist@event_rate)*log(theta),sd=1))
})

