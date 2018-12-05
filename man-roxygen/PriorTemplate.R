#' @param dist data distribution or prior distribution (depending on context) [TODO: make this clearer by renaming arguments in generics!]
#' @param interval numeric vector of length two specifying the interval to condition on.
#' @param x1 stage one outcome (note that n1 is available from \code{design})
#' @param ... further optional arguments
#' @param f univariate function of parameter to get expected value of
#' @param prior prior distribution
#' @param n1 stage-one sample size
