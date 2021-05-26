#' Student's t data distribution
#'
#' Implements exact t-distributions instead of a normal approximation
#'
#' @template DataDistributionTemplate
#'
#' @rdname StudentDataDistribution-class
#' @exportClass Student
setClass("Student", representation(
  two_armed  = "logical",
  multiplier = "numeric"
),
contains = "DataDistribution")


#' @param two_armed logical indicating if a two-armed trial is regarded
#'
#' @examples
#' datadist <- Student(two_armed = TRUE)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname StudentDataDistribution-class
#' @export
Student <- function(two_armed = TRUE) {
  new("Student", two_armed = two_armed, multiplier = ifelse(two_armed, 2, 1))
}


#' @examples
#' probability_density_function(Student(TRUE), 1, 40, 1.1)
#'
#' @rdname probability_density_function
#' @export
setMethod("probability_density_function", signature("Student", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
            return(stats::dt(x, df = max(1, dist@multiplier * (n - 1)), ncp = sqrt(n / dist@multiplier) * theta))
          })


#' @examples
#' cumulative_distribution_function(Student(two_armed = FALSE), .75, 50, .9)
#'
#' @rdname cumulative_distribution_function
#' @export
setMethod("cumulative_distribution_function", signature("Student", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
            return(stats::pt(x, df = max(1, dist@multiplier * (n - 1)), ncp = sqrt(n / dist@multiplier) * theta))
          })



#' @param probs vector of probabilities
#' @rdname StudentDataDistribution-class
#' @export
setMethod("quantile", signature("Student"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
            return(stats::qt(probs, df = max(1, x@multiplier * (n - 1)), ncp = sqrt(n / x@multiplier) * theta))
          })



#' @rdname StudentDataDistribution-class
#'
#' @param object object of class \code{Student}
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("Student", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
            if (!is.null(seed)) set.seed(seed)
            return(stats::rt(nsim, df = max(1, object@multiplier * (n - 1)), ncp = sqrt(n / object@multiplier) * theta))
          })




setMethod("print", signature('Student'), function(x, ...) {
  glue::glue(
    "{class(x)[1]}<{if (x@two_armed) 'two-armed' else 'single-armed'}>"
  )
})
