#' Group-sequential design
#'
#' [ToDo]
#'
#' @slot n1 stage-one sample size
#' @slot c1f early stopping for futility boundary
#' @slot c1e early stopping for efficacy boundary
#' @slot n2 stage-two sample size
#' @slot c2_pivots vector of length order giving the values of c2 at the
#'     pivot points of the numeric integration rule [TODO: these are not available during construction]
#' @slot x1_norm_pivots normalized pivots for integration rule (in [-1, 1])
#' @slot weights weights of conditional score values at x1_norm_pivots for
#'     approximating the integral over x1.
#'
#' @exportClass GSDesign
setClass("GSDesign", representation(
        n1        = "numeric",
        c1f       = "numeric",
        c1e       = "numeric",
        n2        = "numeric",
        c2_pivots = "numeric",
        x1_norm_pivots = "numeric",
        weights   = "numeric"
    ),  contains = "TwoStageDesign"
    )



#' @param n1 cf. slot
#' @param c1f cf. slot
#' @param c1e cf. slot
#' @param n2 cf. slot
#' @param c2_pivots cf. slot
#' @param x1_norm_pivots cf. slot
#' @param weights cf. slot
#' @param ... further optional arguments
#'
#' @rdname GSDesign-class
#' @export
GSDesign <- function(n1, c1f, c1e, n2, c2_pivots, x1_norm_pivots, weights) {
    if (any(diff(sapply(list(c2_pivots, x1_norm_pivots, weights), length)) != 0))
        stop("pivots and weights must all be of the same length")
    if (any(x1_norm_pivots < -1) | any(x1_norm_pivots > 1))
        stop("x1_norm_pivots must be in [-1, 1], is scaled automatically")
    if (any(weights <= 0))
        stop("weights must be positive")
    new("GSDesign", n1 = n1, c1f = c1f, c1e = c1e, n2 = n2,
        c2_pivots = c2_pivots, x1_norm_pivots = x1_norm_pivots, weights = weights)
}



#' @details GQDesign simply creates a TwoStageDesign object with a Gaussian
#'     quadrature numerical integration rule.
#'
#' @param order order (i.e. number of pivot points in the interior of [c1f, c1e])
#'     of the Gaussian quadrature rule to use for integration
#'
#' @rdname GSDesign-class
#' @export
GS_GQDesign <- function(n1, c1f, c1e, n2, c2_pivots, order) {
    if (length(c2_pivots) != order)
        stop("length of pivot vectors does not fit")
    rule <- GaussLegendreRule(order)
    GSDesign(n1 = n1, c1f = c1f, c1e = c1e, n2 = n2, c2_pivots = c2_pivots,
             x1_norm_pivots = rule$nodes, weights = rule$weights)
}



#' @param x object to get parameters from
#'
#' @rdname GSDesign-class
#' @export
setMethod("tunable_parameters", signature("GSDesign"),
          function(x, ...) c(x@n1, x@c1f, x@c1e, x@n2, x@c2_pivots))




#' @param params vector of design parameters (must be in same order as returned
#'     by \code{as.numeric(design)})
#' @param object object to update
#'
#' @rdname GSDesign-class
#' @export
setMethod("update", signature("GSDesign"),
    function(object, params, ...) {
        k <- length(object@weights)
        if( (length(params) - 4) != k)
            stop("parameter length does not fit")
        new("GSDesign",
            n1  = params[1],
            c1f = params[2],
            c1e = params[3],
            n2  = params[4],
            c2_pivots = params[5:(length(params))],
            x1_norm_pivots = object@x1_norm_pivots,
            weights = object@weights)
    })



#' @param x1 stage-one outcome
#' @param d design object
#'
#' @rdname GSDesign-class
#' @export
setMethod("n2", signature("GSDesign", "numeric"),
          function(d, x1, ...) ifelse(x1 < d@c1f | x1 > d@c1e, 0, d@n2) )




#' Convert a group-sequential design to a two-stage design
#'
#' @param d object of class \code{GSDesign}
#'
#' @export

gs2ts <- function(d){
    if(class(d) != "GSDesign")
        stop("d must be of class GSDesign")
    new("TwoStageDesign", n1 = d@n1, c1f = d@c1f, c1e = d@c1e,
                    n2_pivots = rep(d@n2, length(d@weights)),  c2_pivots = d@c2_pivots,
                    x1_norm_pivots = d@x1_norm_pivots, weights = d@weights)
}



#' Find optimal group-sequential design by constraint minimization
#'
#' \code{minimize} takes an unconditioonal score [TODO: need joint superclass for UnconditionalDesign and AffineScores...]
#' and a constraint set (or single constraint) and solves the corresponding
#' constraint minimization problem using \code{nloptr} (using COBYLA by default).
#'
#' @param objective objective function
#' @param subject_to constraint collection
#' @param initial_design initial guess (x0 for nloptr)
#' @param lower_boundary_design design specifying the lower boundary
#' @param upper_boundary_design design specifying the upper boundary
#' @param opts options list passed to nloptr
#' @param ... further optional arguments passed to \code{\link{nloptr}}
#'
#' @export
minimize_gs <- function(objective, subject_to, initial_design,
                        lower_boundary_design, upper_boundary_design,
                        opts = list(
                             algorithm   = "NLOPT_LN_COBYLA",
                            xtol_rel    = 1e-4,
                            maxeval     = 2500
                        ), ...) {

    f_obj <- function(params) evaluate(objective, update(initial_design, params))

    g_cnstr <- function(params) {
        design <- update(initial_design, params)
        user_cnstr <- evaluate(subject_to, design)
        return(c(
            user_cnstr,
            design@c1f - design@c1e + .1, # ensure c1e > c1f
            diff(c2(design, scaled_integration_pivots(design))) # make c2() monotone
        ))
    }

    if (any(g_cnstr(tunable_parameters(initial_design)) > 0))
        warning("initial design is infeasible!")

    res <- nloptr::nloptr(
        x0 = tunable_parameters(initial_design),
        lb = tunable_parameters(lower_boundary_design),
        ub = tunable_parameters(upper_boundary_design),
        eval_f      = f_obj,
        eval_g_ineq = g_cnstr,
        opts = opts,
        ...
    )

    # TODO: error handling

    return(update(initial_design, res$solution))

}
