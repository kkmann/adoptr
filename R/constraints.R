#' Formulating constraints
#'
#' [TODO]
#'
#' [TODO: currently we only support scores on the LHS - this can be easily extended!]
#'
#' [TODO: do we also want to support scores vs scores comparisons?]
#'
#' @param e1 first comparator
#' @param e2 second comparator
#' @template dotdotdotTemplate
#'
#' @exportClass Constraint
setClass("Constraint")

#' @param s constraint to evaluate [TODO make naming of arguments for evaluate more generic!]
#' @param design design to evaluate
#'
#' @rdname Constraint-class
#' @export
setMethod("evaluate", signature("Constraint", "Design"),
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
          function(e1, e2) new("ConditionalConstraint", score = -1 * e1, rhs = - e2))
# TODO: don't be too dogmatic, we can implement all constructors for Scores if
# we define a new abstract class Score...





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
          function(e1, e2) new("UnconditionalConstraint", score = -1 * e1, rhs = - e2))





#' Collection of constraints
#'
#' @slot unconditional_constraints [todo]
#' @slot conditional_constraints [todo]
#'
#' @exportClass ConstraintsCollection
setClass("ConstraintsCollection", representation(
        unconditional_constraints = "list",
        conditional_constraints = "list"))


#' @rdname ConstraintsCollection-class
#' @export
setMethod("evaluate", signature("ConstraintsCollection", "Design"),
          function(s, design, ...) {
              # evaluate the unconditional constraints
              unconditional <- sapply(s@unconditional_constraints, function(cnstr) evaluate(cnstr, design, ...))
              # the conditional constraints must be evaluated 'everywhere'
              # for most designs this will be a set of pivot points in the
              # continuation area
              # TODO: do we want to support conditional constraints
              # on the stopping regions as well? might be handy, e.g. minimal
              # sample size for stopping for efficacy...?
              # to do that we need an additional argument to
              # evaluate(ConditionalScore, Design, where = c(x1, "continuation", "early_futility", "early_efficacy"))
              # which then returns a Design-specific vector evaluating the ConditionalScore on a grid.
              # BETTER: dispatch on class of x1, currently implemented for numeric,
              # for character can implement evaluation on "continuation"
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
    res <- new("ConstraintsCollection", unconditional, conditional)
    return(res)
}
