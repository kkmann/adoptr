#' Formulating constraints
#'
#' Conceptually, constraints work very similar to scores (any score can be put in
#' a constraint).
#' Currently, only constraints of the form 'score <=/>= x' are admissable.#'
#'
#' [TODO: currently we only support scores on the LHS - this can be easily extended!]
#'
#' [TODO: do we also want to support scores vs scores comparisons?]
#'
#' [TODO: conditional constraints are only evaluated at the continuation area]
#'
#' @param e1 first comparator
#' @param e2 second comparator
#' @template dotdotdotTemplate
#'
#' @examples
#' cp          <- ConditionalPower(Normal, PointMassPrior(0.4))
#' pow         <- integrate(cp)
#' constraint1 <- pow >= 0.8 # an unconditional power constraint
#' constraint2 <- cp >= 0.7 # a conditional power constraint
#'
#' @exportClass Constraint
setClass("Constraint")

#' @param s constraint to evaluate
#' @param design TwoStageDesign to evaluate
#'
#' @rdname Constraint-class
#' @export
setMethod("evaluate", signature("Constraint", "TwoStageDesign"),
          function(s, design, ...) evaluate(s@score, design, ...) - s@rhs )





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





#' Collection of constraints
#'
#' @slot unconditional_constraints [todo]
#' @slot conditional_constraints [todo]
#'
#' @param s constraint collection
#' @param design design
#'
#' @exportClass ConstraintsCollection
setClass("ConstraintsCollection", representation(
        unconditional_constraints = "list",
        conditional_constraints = "list"))


#' @rdname ConstraintsCollection-class
#' @export
setMethod("evaluate", signature("ConstraintsCollection", "TwoStageDesign"),
          function(s, design, ...) {
              # TODO: we will want to allow users to chose where the conditional constraints should apply,
              # e.g.  continuation, early efficacy etc.
              x1_cont <- scaled_integration_pivots(design)
              unconditional <- as.numeric(sapply(
                  s@unconditional_constraints,
                  function(cnstr) evaluate(cnstr, design, ...)
              ))
              conditional <- as.numeric(sapply(
                  s@conditional_constraints,
                  function(cnstr) evaluate(cnstr, design, x1_cont, ...)
              ))
              return(c(unconditional, conditional))
        })


#' @description \code{subject_to(...)} can be used to generate an object of class
#'      \code{ConstrintsCollection} from an arbitrary number of (un)conditional
#'      constraints.
#'
#' @param ... arbitrary number of (un)conditional constraints
#'
#' @rdname ConstraintsCollection-class
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
