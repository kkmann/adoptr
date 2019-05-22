#' Formulating Constraints
#'
#' Conceptually, constraints work very similar to scores (any score can be put in
#' a constraint).
#' Currently,  constraints of the form 'score <=/>= x',
#' 'x <=/>= score' and 'score <=/>= score' are admissable.
#'
#' @template s
#' @template design
#' @template optimization
#' @template dotdotdot
#' @template object
#' @param e1 left hand side (score or numeric)
#' @param e2 right hand side (score or numeric)
#'
#' @seealso \code{\link{minimize}}
#'
#' @examples
#' design <- OneStageDesign(50, 1.96)
#'
#' cp     <- ConditionalPower(Normal(), PointMassPrior(0.4, 1))
#' pow    <- Power(Normal(), PointMassPrior(0.4, 1))
#'
#' # unconditional power constraint
#' constraint1 <- pow >= 0.8
#' evaluate(constraint1, design)
#'
#' # conditional power constraint
#' constraint2 <- cp  >= 0.7
#' evaluate(constraint2, design, .5)
#' constraint3 <- 0.7 <= cp # same as constraint2
#' evaluate(constraint3, design, .5)
#'
#' @name Constraints
NULL



setClass("Constraint")
setClass("ConditionalConstraint", representation(
        score = "ConditionalScore",
        rhs   = "numeric"
    ),
    contains = "Constraint")
setClass("UnconditionalConstraint", representation(
        score = "UnconditionalScore",
        rhs   = "numeric"
    ),
    contains = "Constraint")





#' @rdname Constraints
#' @export
setMethod("evaluate", signature("Constraint", "TwoStageDesign"),
          function(s, design, optimization = FALSE, ...) {
              evaluate(s@score, design, optimization, ...) - s@rhs
          })

#' @rdname Constraints
#' @export
setMethod("show", signature(object = "Constraint"),
          function(object) cat(class(object)[1]))





#' @rdname Constraints
#' @export
setMethod("<=", signature("ConditionalScore", "numeric"),
          function(e1, e2) new("ConditionalConstraint", score = e1, rhs = e2))

#' @rdname Constraints
#' @export
setMethod(">=", signature("ConditionalScore", "numeric"),
          function(e1, e2) new("ConditionalConstraint", score = composite({-1*e1}), rhs = -e2))

#' @rdname Constraints
#' @export
setMethod("<=", signature("numeric", "ConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = composite({-1*e2}), rhs = -e1))

#' @rdname Constraints
#' @export
setMethod(">=", signature("numeric", "ConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = e2, rhs = e1))

#' @rdname Constraints
#' @export
setMethod("<=", signature("ConditionalScore", "ConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = composite({e1 - e2}), rhs = 0))

#' @rdname Constraints
#' @export
setMethod(">=", signature("ConditionalScore", "ConditionalScore"),
          function(e1, e2) new("ConditionalConstraint", score = composite({e2 - e1}), rhs = 0))



#' @rdname Constraints
#' @export
setMethod("<=", signature("UnconditionalScore", "numeric"),
          function(e1, e2) new("UnconditionalConstraint", score = e1, rhs = e2))

#' @rdname Constraints
#' @export
setMethod(">=", signature("UnconditionalScore", "numeric"),
          function(e1, e2) new("UnconditionalConstraint", score = composite({-e1}), rhs = -e2))

#' @rdname Constraints
#' @export
setMethod("<=", signature("numeric", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = composite({-e2}), rhs = -e1))

#' @rdname Constraints
#' @export
setMethod(">=", signature("numeric", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = e2, rhs = e1))

#' @rdname Constraints
#' @export
setMethod("<=", signature("UnconditionalScore", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = composite({e1 - e2}), rhs = 0))

#' @rdname Constraints
#' @export
setMethod(">=", signature("UnconditionalScore", "UnconditionalScore"),
          function(e1, e2) new("UnconditionalConstraint", score = composite({e2 - e1}), rhs = 0))




# not user-facing
setClass("ConstraintsCollection", representation(
        unconditional_constraints = "list",
        conditional_constraints   = "list"
    ))



#' Create a collection of constraints
#'
#' \code{subject_to(...)} can be used to generate an object of class
#' \code{ConstraintsCollection} from an arbitrary number of (un)conditional
#' constraints.
#'
#' @param s object of class \code{ConstraintCollection}
#' @template design
#' @template optimization
#' @param ... either constraint objects (for \code{subject_to} or optional arguments passed to \code{evaluate})
#'
#' @return an object of class \code{ConstraintsCollection}
#'
#' @seealso \code{subject_to} is intended to be used for constraint
#'   specification the constraints in \code{\link{minimize}}.
#'
#' @examples
#' # define type one error rate and power
#' toer  <- Power(Normal(), PointMassPrior(0.0, 1))
#' power <- Power(Normal(), PointMassPrior(0.4, 1))
#'
#' # create constrain collection
#' subject_to(
#'   toer  <= 0.025,
#'   power >= 0.9
#' )
#'
#' @aliases ConstraintCollection
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



#' @rdname subject_to
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
