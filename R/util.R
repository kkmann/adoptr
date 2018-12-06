# taken from SMR package - reference properly! Daniel Furtado Ferreira <danielff@dex.ufla.br>
.GaussLegendre <- function(order)
{
    order <- as.integer(order)
    if (order < 0)
        stop("Must be a non-negative number of nodes!")
    if (order == 0)
        return(list(x = numeric(0), w = numeric(0)))
    j   <- 1:(order - 1)
    mu0 <- 2
    b   <- j / (4 * j^2 - 1)^0.5
    A   <- rep(0, order * order)
    A[(order + 1) * (j - 1) + 2] <- b
    A[(order + 1) * j] <- b
    dim(A) <- c(order, order)
    sd <- eigen(A, symmetric = TRUE)
    w  <- rev(as.vector(sd$vectors[1, ]))
    w  <- mu0 * w^2
    x  <- rev(sd$values)
    return(list(
        nodes = x,
        weights = w
    ))
}
