#' Find optimal two-stage design by constraint minimization
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
setGeneric("minimize", function(objective, subject_to, initial_design,
                                lower_boundary_design, upper_boundary_design,
                                opts = list(
                                    algorithm   = "NLOPT_LN_COBYLA",
                                    xtol_rel    = 1e-4,
                                    maxeval     = 2500), ...) standardGeneric("minimize"))
#' @rdname minimize
#' @export
setMethod("minimize", signature("UnconditionalScore", "ConstraintsCollection", "TwoStageDesign"),
     function(objective, subject_to, initial_design,
                         lower_boundary_design, upper_boundary_design,
                         opts, ...) {

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

    })
