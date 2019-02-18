#' Working with scores
#'
#' Both \code{\link{ConditionalScore-class}} as well as
#' \code{\link{UnconditionalScore-class}} implement \code{evaluate} methods
#' which handle the actual computation of the score given a design and
#' (for conditional scores) an interim result.
#' Conditional scores additionally implement an \code{integrate} method
#' to obtain the corresponding unconditional \code{\link{IntegralScore-class}}.
#'
#' @name working-with-scores
NULL
#' @param s score
#' @param design \code{TwoStageDesign} object
#' @template dotdotdotTemplate
#'
#' @rdname working-with-scores
#' @export
setGeneric("evaluate", function(s, design, ...) standardGeneric("evaluate"))


#' @rdname working-with-scores
#' @export
setGeneric("integrate", function(s, ...) standardGeneric("integrate"))


# not user facing
setGeneric(".evaluate", function(s, design, ...) standardGeneric(".evaluate"))



# internal use only
setClass("AbstractConditionalScore")

# needs to implement:
# setMethod("evaluate", signature("AbstractConditionalScore", "TwoStageDesign"),
#          function(s, design, x1, ...) stop("not implemented"))




#' Class for conditional scoring function
#'
#' [ToDo]
#'
#' @template ConditionalScoreTemplate
#'
#' @exportClass ConditionalScore
setClass("ConditionalScore", representation(
        distribution = "DataDistribution",
        prior = "Prior"
    ),
    contains = "AbstractConditionalScore")


#' @describeIn ConditionalScore not implemented, just raises a 'not implemented'
#'     error.
setMethod("evaluate", signature("ConditionalScore", "TwoStageDesign"),
          function(s, design, x1, ...) stop("not implemented"))


#' @describeIn ConditionalScore integrate a \code{ConditionalScore} over the
#'     stage-one outcome; return object of class \code{\link{IntegralScore-class}}.
setMethod("integrate", signature("ConditionalScore"),
          function(s, ...) new("IntegralScore", cs = s) )


#' Score arithmetic
#'
#' To facilitate working with simple weighted sums of scores,
#' \code{otsd} supports some basic arithmetic operations on score object
#' (both conditional and unconditional ones).
#' Scores can be scalar-multiplied by a constant and added to produce new
#' scores.
#' Conditional and unconditional scores cannot be mixed.
#'
#' @param e1 first summand / factor
#' @param e2 second summand / factor
#' @name score-arithmetic
NULL
#' @rdname score-arithmetic
setMethod("+", signature("ConditionalScore", "numeric"),
          function(e1, e2) AffineConditionalScore(list(e1), 1, e2) )
#' @rdname score-arithmetic
setMethod("+", signature("numeric", "ConditionalScore"),
          function(e1, e2) e2 + e1 )
## TODO: check for duplicate scores and combine!
#' @rdname score-arithmetic
setMethod("+", signature("ConditionalScore", "ConditionalScore"),
          function(e1, e2) AffineConditionalScore(list(e1, e2), c(1, 1), 0) )


#' @rdname score-arithmetic
setMethod("*", signature("ConditionalScore", "numeric"),
          function(e1, e2) AffineConditionalScore(list(e1), e2, 0) )
#' @rdname score-arithmetic
setMethod("*", signature("numeric", "ConditionalScore"),
          function(e1, e2) e2 * e1 )





#' Abstract class for unconditional scoring function
#'
#' [ToDo]
#'
#'
#' @exportClass UnconditionalScore
setClass("UnconditionalScore")


#' @param s an \code{IntegralScore}
#' @param design a \code{TwoStageDesign}
#' @template dotdotdotTemplate
#'
#' @rdname UnconditionalScore-class
#' @export
setMethod("evaluate", signature("UnconditionalScore", "TwoStageDesign"),
          function(s, design, ...) stop("not implemented") )


# not user facing
setMethod(".evaluate", signature("UnconditionalScore", "TwoStageDesign"),
          function(s, design, ...) stop("not implemented") )


#' @rdname score-arithmetic
setMethod("+", signature("UnconditionalScore", "numeric"),
          function(e1, e2) AffineUnconditionalScore(list(e1), 1, e2) )

#' @rdname score-arithmetic
setMethod("+", signature("numeric", "UnconditionalScore"),
          function(e1, e2) e2 + e1 )
## TODO: check for duplicate scores and combine!

#' @rdname score-arithmetic
setMethod("+", signature("UnconditionalScore", "UnconditionalScore"),
          function(e1, e2) AffineUnconditionalScore(list(e1, e2), c(1, 1), 0) )


#' @rdname score-arithmetic
setMethod("*", signature("UnconditionalScore", "numeric"),
          function(e1, e2) AffineUnconditionalScore(list(e1), e2, 0) )

#' @rdname score-arithmetic
setMethod("*", signature("numeric", "UnconditionalScore"),
          function(e1, e2) e2 * e1 )




#' Unconditional score class obtained by integration of a \code{ConditionalScore}
#'
#' @param s an \code{IntegralScore}
#' @param design a \code{TwoStageDesign}
#'
#' @slot cs the underlying \code{ConditionalScore}
#'
#' @exportClass IntegralScore
setClass("IntegralScore", representation(
        cs = "ConditionalScore"
    ),
    contains = "UnconditionalScore")


#' @param specific logical, flag for switching to design-specific implementation.
#' @param ... further optimal arguments
#'
#' @describeIn IntegralScore generic implementation of evaluating an integral
#'     score. Uses adaptive Gaussian quadrature for integration and might be
#'     more efficiently implemented by specific \code{TwoStageDesign}-classes
#'     (cf. \code{\link{.evaluate}}).
setMethod("evaluate", signature("IntegralScore", "TwoStageDesign"),
          function(s, design, specific = TRUE, ...) {
              # TODO: currently ignores the possibility of early stopping/uncontinuus
              # conditional scores - might get better when checking for early stopping
              # and integrating separately!
              if (specific) { # use design-specific implementation
                  return(.evaluate(s, design, ...))
              } else {
                  # use generic approach
                  # integrand is the conditional score as function of z1 times the
                  # predictive pdf given the scores prior
                  poef <- predictive_cdf(s@cs@distribution, s@cs@prior, design@c1f, design@n1)
                  poee <- 1 - predictive_cdf(s@cs@distribution, s@cs@prior, design@c1e, design@n1)
                  # continuation region
                  integrand   <- function(x1) evaluate(s@cs, design, x1, ...) *
                      predictive_pdf(s@cs@distribution, s@cs@prior, x1, design@n1, ...)
                  # use adaptive quadrature to integrate - only relies on generic interface
                  # provided by 'TwoStageDesign', no special optimization for particular
                  # design implementation
                  mid_section <- stats::integrate(integrand, design@c1f, design@c1e)$value
                  # compose
                  res <- poef * evaluate( # score is constant on early stopping region
                      s@cs, design,
                      design@c1f - sqrt(.Machine$double.eps) # slightly smaller than stopping for futility
                  ) +
                      mid_section +
                      poee * evaluate(
                          s@cs, design,
                          design@c1e + sqrt(.Machine$double.eps)
                      )
                  return(res)
              }
          })




# not user facing!
setMethod(".evaluate", signature("IntegralScore", "TwoStageDesign"),
          function(s, design, ...) {
              # use design specific implementation tailored to this particular
              # implementation (Gauss Quadrature N points here)
              poef <- predictive_cdf(s@cs@distribution, s@cs@prior, design@c1f, design@n1)
              poee <- 1 - predictive_cdf(s@cs@distribution, s@cs@prior, design@c1e, design@n1)
              # continuation region
              integrand   <- function(x1) evaluate(s@cs, design, x1, ...) *
                  predictive_pdf(s@cs@distribution, s@cs@prior, x1, design@n1, ...)
              mid_section <- integrate_rule(
                  integrand, design@c1f, design@c1e, design@x1_norm_pivots, design@weights
              )
              # compose
              res <- poef * evaluate( # score is constant on early stopping region
                  s@cs, design,
                  design@c1f - sqrt(.Machine$double.eps) # slightly smaller than stopping for futility
              ) +
              mid_section +
              poee * evaluate(
                  s@cs, design,
                  design@c1e + sqrt(.Machine$double.eps)
              )
              return(res)
          })




# not user facing!
setMethod(".evaluate", signature("IntegralScore", "OneStageDesign"),
          function(s, design, ...) {
              # use design specific implementation tailored to this particular
              # implementation (Gauss Quadrature N points here)
              poef <- predictive_cdf(s@cs@distribution, s@cs@prior, design@c1f, design@n1)
              poee <- 1 - predictive_cdf(s@cs@distribution, s@cs@prior, design@c1e, design@n1)
              res  <- poef * evaluate( # score is constant on early stopping region
                      s@cs, design,
                      design@c1f - sqrt(.Machine$double.eps) # slightly smaller than stopping for futility
                  ) +
                  poee * evaluate(
                      s@cs, design,
                      design@c1e + sqrt(.Machine$double.eps)
                  )
              return(res)
          })
