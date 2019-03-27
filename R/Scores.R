#' Evaluation of a score
#'
#' Both \code{\link{ConditionalScore}} as well as
#' \code{\link{UnconditionalScore}} implement \code{evaluate} methods
#' which handle the actual computation of the score given a design and
#' (for conditional scores) an interim result.
#'
#' @seealso Conditional scores additionally implement an \code{expected} method
#'    to obtain the corresponding unconditional \code{\link{IntegralScore}}.
#'
#' @details The method \code{evaluate} is preimplemented for all preimplemented
#'    scores in \pkg{adoptr}.
#'    An example on working with scores is presented
#'    \href{https://kkmann.github.io/adoptr/articles/score-and-constraints-arithmetic.html}{here}.
#'
#' @param s score
#' @param design \code{TwoStageDesign} object
#'
#' @export
setGeneric("evaluate", function(s, design, ...) standardGeneric("evaluate"))


#' Compute the expectation of a conditional score
#'
#' By the method \code{expected} any \code{\link{ConditionalScore}}
#' can be integrated over the full \ifelse{html}{\out{x<sub>1</sub>}}{\eqn{x_1}}-range and returns an
#' \code{\link{IntegralScore}}. I.e., for a conditional score
#' \ifelse{html}{\out{s(design, x<sub>1</sub>)}}{\eqn{s(design, x_1)}}
#' the integral \ifelse{html}{\out{&int; s(design, x<sub>1</sub>) d x<sub>1</sub>}}{\eqn{\int s(design, x_1) d x_1}}
#' is computed.
#'
#' @param s ConditionalScore
#' @template dotdotdot
#'
#' @return an object of class \code{\link{IntegralScore}}
#'
#' @export
setGeneric("expected", function(s, ...) standardGeneric("expected"))


# not user facing
setGeneric(".evaluate", function(s, design, ...) standardGeneric(".evaluate"))



# internal use only
setClass("AbstractConditionalScore")



#' Class for conditional scoring function
#'
#' \code{ConditionalScore} is an abstract class for conditional scores.
#' It requires the two slots \code{distribution} and \code{prior} that
#' determine the data distribution and the prior distribution for the effect
#' parameter. When defining a specific  \code{ConditionalScore}, a corresponding
#' method \code{evaluate} needs to be defined, too.
#' Any \code{ConditionalScore} can be transformed to an unconditional
#' \code{\link{IntegralScore}} by means of the method \code{\link{expected}}.
#'
#' @seealso The common conditional scores \code{\link{ConditionalPower}}
#'    and \code{\link{ConditionalSampleSize}} are preimplemented in \pkg{adoptr}.
#'
#' @examples
#' cp <- ConditionalPower(Normal(), PointMassPrior(0, 1))
#' ep <- expected(cp)
#' design <- TwoStageDesign(50, 0, 2, 50, 2, 5)
#' evaluate(ep, design)
#'
#' @aliases ConditionalScore
#' @exportClass ConditionalScore
setClass("ConditionalScore", representation(
        distribution = "DataDistribution",
        prior = "Prior"
    ),
    contains = "AbstractConditionalScore")


#' @examples
#' # creates power under point mass prior
#' expected(ConditionalPower(Normal(), PointMassPrior(.3, 1)))
#'
#' @rdname expected
#' @export
setMethod("expected", signature("ConditionalScore"),
          function(s, ...) new("IntegralScore", cs = s) )


#' @param object object of class \code{ConditionalScore}
#'
#' @rdname ConditionalScore-class
#' @export
setMethod("show", signature(object = "ConditionalScore"),
          function(object) cat(class(object)[1]))

# Ich würde hier noch ein Beispiel mit minimize() angeben, z.B.
# initial_design <- ...
# opt_res <- minimize(ess, subject_to(power >= 0.8, toer  <= .025), initial_design)
# evaluate(ess, opt_res$design)

#' Score arithmetic
#'
#' To facilitate working with simple weighted sums of scores,
#' \code{\link{adoptr}} supports some basic arithmetic operations on score objects
#' (both conditional and unconditional ones).
#' Scores can be scalar-multiplied by a constant and added to produce new
#' scores.
#' Conditional and unconditional scores cannot be mixed.
#'
#' @examples
#' design <- TwoStageDesign(
#'   n1    = 25,
#'   c1f   = 0,
#'   c1e   = 2.5,
#'   n2    = 50,
#'   c2    = 1.96,
#'   order = 7L
#' )
#' ess <- expected(ConditionalSampleSize(Normal(), PointMassPrior(.4, 1.0)))
#' power <- expected(ConditionalPower(Normal(), PointMassPrior(.4, 1.0)))
#' evaluate(ess + 50*power, design)
#'
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




#' Class for unconditional scoring function
#'
#' \code{UnconditionalScore} is an abstract class for unconditional scores.
#'
#' When defining a new \code{UnconditionalScore}, a corresponding
#' method \code{\link{evaluate}} needs to be defined, too.
#'
#' @examples
#' design <- TwoStageDesign(
#'   n1    = 25,
#'   c1f   = 0,
#'   c1e   = 2.5,
#'   n2    = 50,
#'   c2    = 1.96,
#'   order = 7L
#' )
#'
#' cp <- ConditionalPower(Normal(), PointMassPrior(.1, 1.0))
#' ep <- expected(cp)
#'
#' evaluate(ep, design) # .06081054
#'
#'
#'
#' @seealso There are regularization scores \code{\link{N1}} and
#'    \code{\link{AverageN2}} for sample sizes.
#'    The class \code{\link{IntegralScore}} is a specific subclass that defines
#'    uncondtional scores which are expected conditional scores.
#'
#' @aliases UnconditionalScore
#' @exportClass UnconditionalScore
setClass("UnconditionalScore")



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




#' Unconditional score class obtained by integration of a
#' \code{\link{ConditionalScore}}
#'
#' @slot cs the underlying \code{\link{ConditionalScore}}
#'
#' @examples
#' expected(ConditionalPower(Normal(), PointMassPrior(.4, 1.0)))
#'
#' @seealso The method \code{\link{expected}} creates a \code{IntegralScore}
#'    from a \code{\link{ConditionalScore}}.
#'
#' @aliases IntegralScore
#' @exportClass IntegralScore
setClass("IntegralScore", representation(
    cs = "ConditionalScore"
),
contains = "UnconditionalScore")


#' @param object object of class \code{IntegralScore}
#'
#' @rdname IntegralScore-class
#' @export
setMethod("show", signature(object = "IntegralScore"),
          function(object) cat(class(object)[1]))


#' @examples
#' # create a dummy design
#' design <- TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L)
#'
#' # define type one error als IntegralScore
#' toer <- expected(ConditionalPower(Normal(), PointMassPrior(.0, 1)))
#'
#' # evaluate
#' evaluate(toer, design)
#'
#' @template optimization
#' @param subdivisions integer, number of subdivisions that is used for integration
#'    on the continuation region; results become more precise with increased
#'    number of subdivisions; default is \code{10000L}.
#' @template dotdotdot
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("IntegralScore", "TwoStageDesign"),
          function(s, design, optimization = FALSE, subdivisions = 10000L, ...) {
              if (optimization) { # use design-specific implementation
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
                  mid_section <- stats::integrate(
                      integrand, design@c1f, design@c1e, subdivisions = subdivisions)$value
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
              # probability of early futility
              poef <- predictive_cdf(
                  s@cs@distribution, s@cs@prior, design@c1f, n1(design, round = FALSE))
              # probability of early efficacy
              poee <- 1 - predictive_cdf(
                  s@cs@distribution, s@cs@prior, design@c1e, n1(design, round = FALSE))
              # integrand: conditional score times predictive PDF
              integrand <- function(x1) {
                  evaluate(s@cs, design, x1, optimization = TRUE, ...) *
                      predictive_pdf(s@cs@distribution, s@cs@prior, x1, n1(design, round = FALSE), ...)
              }
              mid_section <- integrate_rule(
                  integrand,
                  design@c1f, design@c1e,
                  design@x1_norm_pivots,
                  design@weights
              )
              # compose
              res <- poef * evaluate( # score is constant on early stopping region
                  s@cs, design,
                  design@c1f - sqrt(.Machine$double.eps), # slightly smaller than stopping for futility
                  optimization = TRUE,
                  ...
              ) +
                  mid_section +
                  poee * evaluate(
                      s@cs, design,
                      design@c1e + sqrt(.Machine$double.eps),
                      optimization = TRUE,
                      ...
                  )
              return(res)
          })




# not user facing!
setMethod(".evaluate", signature("IntegralScore", "OneStageDesign"),
          function(s, design, ...) {
              # use design specific implementation tailored to this particular
              # implementation (Gauss Quadrature N points here)
              poef <- predictive_cdf(
                  s@cs@distribution, s@cs@prior, design@c1f, n1(design, round = FALSE))
              poee <- 1 - predictive_cdf(
                  s@cs@distribution, s@cs@prior, design@c1e, n1(design, round = FALSE))
              res  <- poef * evaluate( # score is constant on early stopping region
                  s@cs, design,
                  design@c1f - sqrt(.Machine$double.eps), # slightly smaller than stopping for futility
                  optimization = TRUE,
                  ...
              ) +
                  poee * evaluate(
                      s@cs, design,
                      design@c1e + sqrt(.Machine$double.eps),
                      optimization = TRUE,
                      ...
                  )
              return(res)
          })



