#' Formulating constraints
#'
#' Conceptually, constraints work very similar to scores (any score can be put in
#' a constraint).
#' Currently,  constraints of the form 'score <=/>= x',
#' 'x <=/>= score' and 'score <=/>= score' are admissable.
#'
#' @param e1 first comparator
#' @param e2 second comparator
#' @template dotdotdot
#'
#' @examples
#' cp          <- ConditionalPower(Normal(), PointMassPrior(0.4, 1))
#' pow         <- expected(cp)
#' constraint1 <- pow >= 0.8 # an unconditional power constraint
#' constraint2 <- cp >= 0.7 # a conditional power constraint
#' constraint3 <- 0.7 <= cp # yields the same as constraint2
#'
#' @aliases Constraint
#' @exportClass Constraint
setClass("Constraint")

#' @examples
#' evaluate(
#'    expected(ConditionalPower(Normal(), PointMassPrior(.0, 1))) <= 0.05,
#'    OneStageDesign(50, 1.96)
#' ) # -0.025
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("Constraint", "TwoStageDesign"),
          function(s, design, optimization = FALSE, ...) {
              evaluate(s@score, design, optimization, ...) - s@rhs
          })



#' @rdname Constraint-class
#'
#' @param object object of class \code{Constraint}
#' @export
setMethod("show", signature(object = "Constraint"),
          function(object) cat(class(object)[1]))



#' @rdname Constraint-class
#' @exportClass ConditionalConstraint
setClass("ConditionalConstraint", representation(
        score = "AbstractConditionalScore",
        rhs   = "numeric"
    ),
    contains = "Constraint")


#' @rdname Constraint-class
#' @export
setMethod("<=", signature("AbstractConditionalScore", "numeric"),
          function(e1, e2) new("ConditionalConstraint", score = e1, rhs = e2))
#' @rdname Constraint-class
#' @export
setMethod(">=", signature("AbstractConditionalScore", "numeric"),
          function(e1, e2) new("ConditionalConstraint", score = -1 * e1, rhs = -e2))
#' @rdname Constraint-class
#' @export
setMethod("<=", signature("numeric", "AbstractConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = -1 * e2, rhs = -e1))
#' @rdname Constraint-class
#' @export
setMethod(">=", signature("numeric", "AbstractConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = e2, rhs = e1))
#' @rdname Constraint-class
#' @export
setMethod("<=", signature("AbstractConditionalScore", "AbstractConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = e1 + (-1) * e2, rhs = 0))
#' @rdname Constraint-class
#' @export
setMethod(">=", signature("AbstractConditionalScore", "AbstractConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = e2 + (-1) * e1, rhs = 0))





#' @rdname Constraint-class
#' @exportClass UnconditionalConstraint
setClass("UnconditionalConstraint", representation(
        score = "UnconditionalScore",
        rhs   = "numeric"
    ),
    contains = "Constraint")


#' @rdname Constraint-class
#' @export
setMethod("<=", signature("UnconditionalScore", "numeric"),
          function(e1, e2) new("UnconditionalConstraint", score = e1, rhs = e2))
#' @rdname Constraint-class
#' @export
setMethod(">=", signature("UnconditionalScore", "numeric"),
          function(e1, e2) new("UnconditionalConstraint", score = -1 * e1, rhs = -e2))
#' @rdname Constraint-class
#' @export
setMethod("<=", signature("numeric", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = -1 * e2, rhs = -e1))
#' @rdname Constraint-class
#' @export
setMethod(">=", signature("numeric", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = e2, rhs = e1))
#' @rdname Constraint-class
#' @export
setMethod("<=", signature("UnconditionalScore", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = e1 + (-1) * e2, rhs = 0))
#' @rdname Constraint-class
#' @export
setMethod(">=", signature("UnconditionalScore", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = e2 + (-1) * e1, rhs = 0))





#' Collection of constraints
#'
#' @slot unconditional_constraints a list of elements of class \code{UnconditionalConstraint}
#' @slot conditional_constraints a list of elements of class \code{ConditionalConstraint}
#'
#' A \code{ConstraintsCollection} is a collection of unconditional and
#' conditional constraints. In order to evaluate these correctly, they
#' have to be defined in two different slots.
#'
#' @seealso A \code{ConstraintsCollection} can be created by \code{\link{subject_to}}.
#'
#' @param s constraint collection
#' @param design design
#'
#' @examples
#' # Define Type one error rate
#' toer <- expected(ConditionalPower(Normal(), PointMassPrior(0.0, 1)))
#'
#' # Define Power at delta = 0.4
#' pow <- expected(ConditionalPower(Normal(), PointMassPrior(0.4, 1)))
#'
#' # Define expected sample size at delta = 0.4
#' ess <- expected(ConditionalSampleSize(Normal(), PointMassPrior(0.4, 1)))
#'
#' # Compute design minimizing ess subject to power and toer constraints
#' \dontrun{
#' minimize(
#'    ess,
#'    subject_to(
#'       toer <= 0.025,
#'       pow  >= 0.9
#'    ),
#'    initial_design = TwoStageDesign(50, .0, 2.0, 60.0, 2.0, 5L)
#' )
#' }
#'
#'
#'
#'
#' @aliases ConstraintsCollection
#' @exportClass ConstraintsCollection
setClass("ConstraintsCollection", representation(
        unconditional_constraints = "list",
        conditional_constraints = "list"))


#' @examples
#' # define power at delta = 0.3 and type one error rate
#' pow  <- expected(ConditionalPower(Normal(), PointMassPrior(.3, 1)))
#' toer <- expected(ConditionalPower(Normal(), PointMassPrior(.0, 1)))
#' # evaluate if power >= 0.8 and toer <= 0.025
#' evaluate(
#'    subject_to(
#'       pow  >= 0.8,
#'       toer <= 0.025
#'    ),
#'    TwoStageDesign(50.0, 0.0, 2.0, rep(60.0, 5), seq(2.0, 0.0, length.out = 5))
#' )
#'
#'
#' @rdname evaluate
#' @export
setMethod("evaluate", signature("ConstraintsCollection", "TwoStageDesign"),
          function(s, design, optimization = FALSE, ...) {
              x1_cont <- scaled_integration_pivots(design)
              unconditional <- as.numeric(sapply(
                  s@unconditional_constraints,
                  function(cnstr) evaluate(cnstr, design, optimization, ...)
              ))
              conditional <- as.numeric(sapply(
                  s@conditional_constraints,
                  function(cnstr) evaluate(cnstr, design, x1_cont, optimization, ...)
              ))
              return(c(unconditional, conditional))
        })


#' Create a collection of constraints
#'
#' \code{subject_to(...)} can be used to generate an object of class
#' \code{\link{ConstraintsCollection}} from an arbitrary number of (un)conditional
#' constraints.
#'
#' @param ... arbitrary number of (un)conditional constraints
#'
#' @return an object of class \code{\link{ConstraintsCollection}}
#'
#' @seealso \code{subject_to} can be used for the constraints in
#'    \code{\link{minimize}}.
#'
#' @examples
#' # Define Type one error rate
#' toer <- expected(ConditionalPower(Normal(), PointMassPrior(0.0, 1)))
#'
#' # Define Power at delta = 0.4
#' pow <- expected(ConditionalPower(Normal(), PointMassPrior(0.4, 1)))
#'
#' # Define expected sample size at delta = 0.4
#' ess <- expected(ConditionalSampleSize(Normal(), PointMassPrior(0.4, 1)))
#'
#' # Compute design minimizing ess subject to power and toer constraints
#' \dontrun{
#' minimize(
#'    ess,
#'    subject_to(
#'       toer <= 0.025,
#'       pow  >= 0.9
#'    ),
#'    initial_design = TwoStageDesign(50, .0, 2.0, 60.0, 2.0, 5L)
#' )
#' }
#' @export
subject_to <- function(...) {
    args <- list(...)
    # sort arguments to conditional vs. unconditional
    conditional <- list()
    unconditional <- list()
    for (i in 1:length(args)) {
        if (is(args[[i]], "ConditionalConstraint")) {
            conditional <- append(conditional, args[i])
        } else {
            if (is(args[[i]], "UnconditionalConstraint")) {
                unconditional <- append(unconditional, args[i])
            } else {
                stop("arguments must be of class ConditionalConstraint or UnconditionalConstraint")
            }
        }
    }
    res <- new("ConstraintsCollection", unconditional_constraints = unconditional, conditional_constraints = conditional)
    return(res)
}
