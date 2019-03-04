#' Compute weights an knots for Gauss-Legendre integration
#'
#' Implementation of Gauss Legendre integration rule and computation of integral
#' based on implementation at https://rdrr.io/cran/SMR/src/R/GaussLegendre.R
#' \code{GaussLegendreRule} gives the nodes and weights for the interval \eqn{(-1, 1)}.
#'
#' @param order Number of nodes to be used
#'
#' @rdname GaussLegendreRule
#' @keywords internal
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


#' Integration via the Gauss-Legendre quadrature
#'
#' @param f function to be integrated
#' @param low lower bound of integration
#' @param up upper bound of integration
#' @param x normalized integration knots on [-1, 1]
#' @param weights integration weights
#'
#' @rdname GaussLegendreIntegral
#' @keywords internal
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




#' Boundary designs
#'
#' Compute boundary designs. These can be used for optimization and are
#' implemented by default in \link{minimize}.
#'
#' @name boundary-designs
#'
#' @param objective An object of class \code{\link{UnconditionalScore-class}} to be minimized
#' @param subject_to A \code{\link{ConstraintsCollection-class}} defining the constraints
#' @param initial_design The initial design
#'
#' @rdname boundary-designs
#' @export
setGeneric("get_lower_boundary_design",
           function(objective, subject_to, initial_design, ...)  standardGeneric("get_lower_boundary_design"))

#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "OneStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              OneStageDesign(5, 0.0)
})


#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "GroupSequentialDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

                  os_design <- minimize(
                      objective      = objective,
                      subject_to     = cnstrs,
                      initial_design = OneStageDesign(200, 2.0)
                      )$design


                      lb_design <- GroupSequentialDesign(5,
                                                         -1,
                                                         os_design@c1f - .05,
                                                         2,
                                                         rep(-2, length(initial_design@c2_pivots)))

                return(lb_design)
})



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "TwoStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

              os_design <- minimize(
                  objective      = objective,
                  subject_to     = cnstrs,
                  initial_design = OneStageDesign(200, 2.0)
              )$design


              lb_design <- TwoStageDesign(5,
                                          -2,
                                          os_design@c1f - .05,
                                          rep(2, length(initial_design@c2_pivots)),
                                          rep(-2, length(initial_design@c2_pivots)))

              return(lb_design)
})




#' @rdname boundary-designs
#' @export
setGeneric("get_upper_boundary_design",
           function(objective, subject_to, initial_design, ...)  standardGeneric("get_upper_boundary_design"))


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "OneStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              OneStageDesign(10000, 5.0)
})


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "GroupSequentialDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

              os_design <- minimize(
                  objective      = objective,
                  subject_to     = cnstrs,
                  initial_design = OneStageDesign(200, 2.0)
              )$design


              ub_design <- GroupSequentialDesign(os_design@n1,
                                                 os_design@c1f + .05,
                                                 5,
                                                 2*os_design@n1,
                                                 rep(5, length(initial_design@c2_pivots)))

              return(ub_design)
})


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design",
          signature("UnconditionalScore", "ConstraintsCollection", "TwoStageDesign"),
          function(objective, subject_to, initial_design, ...) {
              cnstrs <- subject_to
              cnstrs@conditional_constraints <- list()

              os_design <- minimize(
                  objective      = objective,
                  subject_to     = cnstrs,
                  initial_design = OneStageDesign(200, 2.0)
              )$design


              ub_design <- TwoStageDesign(os_design@n1,
                                          os_design@c1f + .05,
                                          5,
                                          rep(2*os_design@n1, length(initial_design@c2_pivots)),
                                          rep(5, length(initial_design@c2_pivots)))

              return(ub_design)
})
