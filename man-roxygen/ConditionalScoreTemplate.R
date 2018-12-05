#' @slot distribution \code{DataDistribution} object specifying the data distribution
#'     given the parameter \eqn{\theta} of the model
#' @slot prior \code{Prior} object specifying the data distribution
#'     given the parameter \eqn{\theta} of the model
#'
#' @param s conditional score object to evaluate
#' @param design a \code{Design}
#' @param x1 stage one outcome (note that n1 is available from \code{design})
#' @param ... further optional parameters
