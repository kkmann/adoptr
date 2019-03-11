#' Regularization via L1 norm
#'
#' Implements the L1-norm of the design's stage-two sample size function.
#' The average of the stage-two sample size without weighting with
#' the data distribution is computed.
#' This can be interpreted as integration over a unifrom prior on
#' the continuation region.
#' Adding the L1-norm with a small regularization parameter can make
#' the optimization more stable.
#'
#' @seealso \code{\link{N1}} for penalizing n1 values
#'
#' @aliases AverageN2
#' @exportClass AverageN2
setClass("AverageN2", representation(
    dummy = 'logical'
    ),
    contains = "UnconditionalScore")


#' @examples
#' avn2 <- AverageN2()
#'
#' @return an object of class \code{\link{AverageN2}}
#'
#' @rdname AverageN2-class
#' @export
AverageN2 <- function() new("AverageN2", dummy = FALSE)


#' @examples
#' evaluate(
#'    AverageN2(),
#'    TwoStageDesign(100, 0.5, 1.5, 60.0, 1.96, order = 5L)
#' ) # 60
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("AverageN2", "TwoStageDesign"),
          function(s, design, optimization = FALSE, subdivisions = 10000L, ...) {
              if (optimization) {
                  # use design-specific implementation
                  return(.evaluate(s, design, ...))
              } else {
                  # generic integration
                  res <- stats::integrate(
                      function(x) n2(design, x, round = TRUE),
                      design@c1f,
                      design@c1e,
                      subdivisions = subdivisions,
                      ...
                  )$value
                  res <- res / (design@c1e - design@c1f)
                  return(res)
              }
          }
)


# not user facing!
setMethod(".evaluate", signature("AverageN2", "TwoStageDesign"),
          function(s, design, ...) {
              # use non-rounded version
              integrate_rule(
                  function(x) n2(design, x, round = FALSE, ...),
                  design@c1f,
                  design@c1e,
                  design@x1_norm_pivots,
                  design@weights
                ) / (design@c1e - design@c1f)
          }
)



#' Quadratic smoothness penalty term
#'
#' \code{SmoothnessN2} is a generic class for implementing a smoothness penalty
#' via the average squared second derivative of the stage two sample size
#' function \eqn{n_2}{n2} of a two-stage design.
#' The only parameter is the width used for the finite differences, \eqn{h}{h}.
#' The generic implementation only evaluates \eqn{n_2}{n2} in the interior of
#' the continuation region of a design.
#'
#' @slot h positive number giving the width of the central finite difference
#'     interval for approximating the second derivative.
#'
#' @aliases SmothnessN2
#' @exportClass SmoothnessN2
setClass("SmoothnessN2", representation(
    h = "numeric"
    ),
    contains = "UnconditionalScore")



#' @param h positive number, see slot \code{h}
#'
#' @return an object of class \code{\link{SmoothnessN2}}
#'
#' @rdname SmoothnessN2-class
#' @export
SmoothnessN2 <- function(h = .1)
    new("SmoothnessN2", h = h)




#' @examples
#' evaluate(
#'    SmoothnessN2(),
#'    TwoStageDesign(50, 0, 2, rep(50, 7), rep(2, 7))
#' ) # 0
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("SmoothnessN2", "TwoStageDesign"),
          function(s, design, optimization = FALSE, subdivisions = 10000L, ...) {
              if (optimization) {
                  # use design-specific implementation
                  return(.evaluate(s, design, ...))
              } else {
                  # use generic approach
                  # integrand is the finite difference approximation of the
                  # squared second derivative
                  integrand <- function(x1) {
                      ((n2(design, x1 + s@h, round = TRUE) - 2 * n2(design, x1, round = TRUE) + n2(design, x1 - s@h, round = TRUE)) / s@h^2)^2
                  }
                  x1_bounds <- c(design@c1f + s@h, design@c1e - s@h)
                  # use adaptive quadrature to integrate - only relies on generic interface
                  # provided by 'Design', no special optimization for particular
                  # design implementation
                  return(1 / diff(x1_bounds) *
                             stats::integrate(integrand, x1_bounds[1], x1_bounds[2], subdivisions = subdivisions, ...)$value)
              }
          }
)

# specific method for class TwoStageDesign
setMethod(".evaluate", signature("SmoothnessN2", "TwoStageDesign"),
          function(s, design, ...){
              leng <- length(design@n2_pivots) # length of vector
              dif1 <- diff(design@n2_pivots[-leng]) # increments of n2
              dif2 <- diff(design@n2_pivots[-1]) # increments of n2
              piv  <- diff(scaled_integration_pivots(design)) # increments of pivots
              # Approximate L2 norm of second derivative
              res <- mean(((-dif1 + dif2) / (piv[-leng + 1] * piv[-1]))^2)
              return(res)
          }
)


#' @examples
#' evaluate(
#'    SmoothnessN2(),
#'    GroupSequentialDesign(50, 0, 2, 50, rep(2, 7))
#' ) # 0
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("SmoothnessN2", "GroupSequentialDesign"),
          function(s, design, ...) 0 )



#' @rdname evaluate
#' @export
setMethod("evaluate", signature("SmoothnessN2", "OneStageDesign"),
          function(s, design, ...) 0 )




#' Regularize n1
#'
#' \code{N1} is a class that computes the \code{n1} value of a design.
#'
#' @aliases N1
#' @exportClass N1
setClass("N1", representation(
    dummy = 'logical'
),
contains = "UnconditionalScore")

#' @examples
#' n1_score <- N1()
#'
#' @return an object of class \code{\link{N1}}
#'
#' @rdname N1-class
#' @export
N1 <- function() new("N1", dummy = FALSE)


#' @examples
#' evaluate(
#'    N1(),
#'    TwoStageDesign(70, 0, 2, rep(60, 6), rep(1.7, 6))
#' ) # 70
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("N1", "TwoStageDesign"),
          function(s, design, optimization = FALSE, ...)
              n1(design, round = !optimization)
          )

