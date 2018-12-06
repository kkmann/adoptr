#' Implementation of Gauss Legendre integration rule and computation of integral
#' based on implementation at https://rdrr.io/cran/SMR/src/R/GaussLegendre.R
#' \code{GaussLegendreRule} gives the nodes and weights for the interval \eqn{(-1, 1)}.
#'
#'
#' @param order Number of nodes to be used
#'
#' @rdname GaussLegendreRule
#' @export

GaussLegendreRule <- function(order) {
    order <- as.integer(order)
    if (order < 0)
        stop("Must be a non-negative number of nodes!")
    if (order == 0)
        return(list(x = numeric(0), w = numeric(0)))
    i  <- 1:order
    j   <- 1:(order-1)
    mu0 <- 2
    b <- j / (4 * j^2 - 1)^0.5
    A <- rep(0, order * order)
    A[(order + 1) * (j - 1) + 2] <- b
    A[(order + 1) * j] <- b
    dim(A) <- c(order, order)
    sd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(sd$vectors[1, ]))
    w <- mu0 * w^2
    x <- rev(sd$values)
    return(data.frame(nodes = x, weights = w))
}


#' Integration via the Gauss Legendre quadrature
#'
#' @param f function to be integrated
#' @param low lower bound of integration
#' @param up upper bound of integration
#' @param rule integration rule as computed by \code{GaussLegendreRule}
#'
#' @rdname GaussLegendreIntegral
#' @export

GaussLegendreIntegral <- function(f, low, up, rule) {
  if(ncol(rule) != 2)
    stop("rule needs to contain two columns: nodes and weights!")
  a <- (up - low) / 2
  b <- (up + low) / 2
  ff <- sapply(rule$nodes, function(x) f(a * x + b))
  res <- a * sum(rule$w * ff)
  return(res)
}
