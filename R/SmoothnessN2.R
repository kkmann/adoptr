#' Quadratic smoothness penalty term
#'
#' \code{SmoothnessN2} is a generic class for implementing a smoothness penalty
#' via the average squared second derivative of the stage two sample size function
#' \code{n2} of a two-stage design.
#' The only parameter is the width used for the finite differences, \code{h}.
#' The generic implementation only evluates \code{n2} in the interior of the
#' continuation region of a design.
#'
#' @slot distribution Data distribution
#' @slot prior Prior distribution
#' @slot h positive number giving the width of the central finite difference
#'     interval for approximating the second derivative.
#'
#' @exportClass SmoothnessN2
setClass("SmoothnessN2", representation(
        distribution = "DataDistribution",
        prior = "Prior",
        h = "numeric"
    ),
    contains = "ConditionalScore")



#' @param distribution see slot
#' @param prior see slot
#' @param h positive number, see slot \code{h}
#'
#' @rdname SmoothnessN2-class
#' @export
SmoothnessN2 <- function(distribution,
                         prior = ContinuousPrior(
                             pdf = function(x) stats::dunif(x, -5, 5),
                             support = c(-5, 5)),
                         h = .1)
    new("SmoothnessN2", distribution = distribution, prior = prior, h = h)




#' A generic implementation for arbitrary two-stage designs based on adaptive
#' Gaussian quadrature integration of the finite-differences approximation to
#' the second derivative is provided.
#' Custom subclasses of \code{Design} might implement this slightly different
#' (cf. [TODO: link to GQDesign implementation]).
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
                  integrand <- function(x1) ((n2(design, x1 + s@h) - 2 * n2(design, x1) + n2(design, x1 - s@h)) / s@h^2)^2
                  x1_bounds <- c(design@c1f + s@h, design@c1e - s@h)
                  # use adaptive quadrature to integrate - only relies on generic interface
                  # provided by 'Design', no special optimization for particular
                  # design implementation
                  return(1 / diff(x1_bounds) *
                             stats::integrate(integrand, x1_bounds[1], x1_bounds[2], subdivisions = 1000)$value)
              }
          })

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
              } )



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

