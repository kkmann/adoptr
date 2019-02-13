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
#' @param s a score of class \code{AverageN2}
#' @param design a \code{TwoStageDesign}
#'
#' @rdname AverageN2-class
#'
#' @exportClass AverageN2
setClass("AverageN2", representation(
    dummy = 'logical'
    ),
    contains = "UnconditionalScore")



#' @rdname AverageN2-class
#' @export
AverageN2 <- function() new("AverageN2", dummy = FALSE)

#' Evaluation of a AverageN2 score
#'
#' @param specific logical, flag for switching to design-specific implementation.
#' @param ... further optimal arguments
#'
#' @describeIn AverageN2 generic implementation of evaluating a smoothness
#'     score. Uses adaptive Gaussian quadrature for integration and might be
#'     more efficiently implemented by specific \code{TwoStageDesign}-classes
#'     (cf. \code{\link{.evaluate}}).
setMethod("evaluate", signature("AverageN2", "TwoStageDesign"),
          function(s, design, specific = TRUE, ...) {
              if (specific) { # use design-specific implementation
                  return(.evaluate(s, design, ...))
              } else {
                  res <- stats::integrate(
                      function(x) n2(design, x),
                      design@c1f,
                      design@c1e
                  )$value
                  res <- res / (design@c1e - design@c1f)
                  return(res)
              }
          }
)


# not user facing!
setMethod(".evaluate", signature("AverageN2", "TwoStageDesign"),
          function(s, design, ...) {
              integrate_rule(
                  function(x) n2(design, x),
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
#' via the average squared second derivative of the stage two sample size function
#' \code{n2} of a two-stage design.
#' The only parameter is the width used for the finite differences, \code{h}.
#' The generic implementation only evluates \code{n2} in the interior of the
#' continuation region of a design.
#'
#' @slot h positive number giving the width of the central finite difference
#'     interval for approximating the second derivative.
#'
#' @exportClass SmoothnessN2
setClass("SmoothnessN2", representation(
    h = "numeric"
    ),
    contains = "UnconditionalScore")



#' @param h positive number, see slot \code{h}
#'
#' @rdname SmoothnessN2-class
#' @export
SmoothnessN2 <- function(h = .1)
    new("SmoothnessN2", h = h)




#' A generic implementation for arbitrary two-stage designs based on adaptive
#' Gaussian quadrature integration of the finite-differences approximation to
#' the second derivative is provided.
#' For \code{TwoStageDesign} and its subclasses a specific implementation
#' is avaible which is recommended to use.
#'
#' @param s an object of class \code{SmoothnessN2}
#' @param design the design to compute the smoothness term for
#' @param specific should a specific implementation be used?
#' @template dotdotdotTemplate
#'
#' @rdname SmoothnessN2-class
#' @export
setMethod("evaluate", signature("SmoothnessN2", "TwoStageDesign"),
          function(s, design, specific = TRUE, ...) {
              if (specific) {
                  # use design-specific implementation
                  return(.evaluate(s, design, ...))
              } else {
                  # use generic approach
                  # integrand is the finite difference approximation of the
                  # squared second derivative
                  integrand <- function(x1) {
                      ((n2(design, x1 + s@h) - 2 * n2(design, x1) + n2(design, x1 - s@h)) / s@h^2)^2
                  }
                  x1_bounds <- c(design@c1f + s@h, design@c1e - s@h)
                  # use adaptive quadrature to integrate - only relies on generic interface
                  # provided by 'Design', no special optimization for particular
                  # design implementation
                  return(1 / diff(x1_bounds) *
                             stats::integrate(integrand, x1_bounds[1], x1_bounds[2])$value)
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
              res <- mean(((-dif1 + dif2) / (piv[-leng+1] * piv[-1]))^2)
              return(res)
          }
)



#' Return smootheness of a group-sequential design as 0.
#'
#' @param s an object of class \code{SmoothnessN2}
#' @param design an object of class \code{GSDesign}
#'
#' @rdname GSDesign-class
#' @export
setMethod("evaluate", signature("SmoothnessN2", "GSDesign"),
          function(s, design, ...) 0 )



#' Return smootheness of a one-stage design as 0.
#'
#' @param s an object of class \code{SmoothnessN2}
#' @param design an object of class \code{OneStageDesign}
#'
#' @rdname OneStageDesign-class
#' @export
setMethod("evaluate", signature("SmoothnessN2", "OneStageDesign"),
          function(s, design, ...) 0 )




#' Penalize n1
#'
#' \code{PenaltyN1} is a class that penalizes the \code{n1} value of
#' a design.
#'
#'
#' @exportClass PenaltyN1
setClass("PenaltyN1", representation(
    dummy = 'logical'
),
contains = "UnconditionalScore")



#' @rdname PenaltyN1-class
#' @export
PenaltyN1 <- function() new("PenaltyN1", dummy = FALSE)


#' Returns the n1-value of a \code{TwoStageDesign}.
#' This can be used as penalty term in the optimization.
#'
#' @param s an object of class \code{PenaltyN1}
#' @param design the design to compute the smoothness term for
#' @template dotdotdotTemplate
#'
#' @rdname PenaltyN1-class
#' @export
setMethod("evaluate", signature("PenaltyN1", "TwoStageDesign"),
          function(s, design, ...) {
              design@n1
          }
)

