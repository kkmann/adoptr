#' @param order \code{integer}, integration order of the employed Gaussian quadrature
#' integration rule to evaluate scores. Automatically set to \code{length(n2_pivots)} if \cr
#' \code{length(n2_pivots) == length(c2_pivots) > 1}, otherwise c2 and n2
#' are taken to be constant in stage-two and replicated to match the number of
#' pivots specified by \code{order}
