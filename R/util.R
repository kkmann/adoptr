# Compute weights an knots for Gauss-Legendre integration
#
# Implementation of Gauss Legendre integration rule and computation of integral
# based on implementation at https://rdrr.io/cran/SMR/src/R/GaussLegendre.R
# \code{GaussLegendreRule} gives the nodes and weights for the interval \eqn{(-1, 1)}.
# internal
GaussLegendreRule <- function(order) {
    order <- as.integer(order)
    if (order < 2)
        stop("At least two nodes are necessary for integration!")
    j   <- 1:(order - 1)
    mu0 <- 2
    b   <- j / (4 * j^2 - 1)^0.5
    A   <- rep(0, order * order)
    A[(order + 1) * (j - 1) + 2] <- b
    A[(order + 1) * j] <- b
    dim(A) <- c(order, order)
    sd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(sd$vectors[1, ]))
    w <- mu0 * w^2
    x <- rev(sd$values)
    return(data.frame(nodes = x, weights = w))
}


# integration via the Gauss-Legendre quadrature, internal
integrate_rule <- function(f, low, up, x, weights) {

  if (!(length(weights) == length(x)))
      stop("x and weights must be of same length")

  if (any(x < -1) | any(x > 1))
    stop("x must be in [-1, 1], is scaled automatically")

  if (any(weights <= 0))
    stop("weights must be positive")

  a  <- (up - low) / 2
  b  <- (up + low) / 2
  ff <- sapply(x, function(x) f(a * x + b))

  return(a * sum(weights * ff))

}
