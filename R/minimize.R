#' Find optimal two-stage design by constraint minimization
#'
#' \code{minimize} takes an unconditioonal score
#' and a constraint set (or single constraint) of conditional and/or unconditional
#' scores and solves the corresponding constraint minimization problem
#' using \code{nloptr} (using COBYLA by default).
#' An initial design has to be defined. It is also possible to defined
#' lower- and upper-boundary designs. If this is not done, these
#' are computed automatically.
#'
#' @param objective objective function
#' @param subject_to constraint collection
#' @param initial_design initial guess (x0 for nloptr)
#' @param lower_boundary_design design specifying the lower boundary.
#' @param upper_boundary_design design specifying the upper boundary
#' @param opts options list passed to nloptr
#' @param ... further optional arguments passed to \code{\link{nloptr}}
#'
#' @return \item{design}{ The resulting optimal design}
#'         \item{nloptr_return}{ Output of the corresponding nloptr call}
#'         \item{call_args}{ The arguments given to the optimization call}
#'
#' @export
minimize <- function(
    objective,
    subject_to,
    initial_design,
    lower_boundary_design = get_lower_boundary_design(objective, subject_to, initial_design),
    upper_boundary_design = get_upper_boundary_design(objective, subject_to, initial_design),
    opts         =  list(
        algorithm   = "NLOPT_LN_COBYLA",
        xtol_rel    = 1e-5,
        maxeval     = 10000
    ),
    ...
) {

    args <- c(as.list(environment()), list(...))

    f_obj <- function(params) {
        evaluate(
            objective,
            update(initial_design, params),
            optimization = TRUE # evaluate in optimization context!
        )
    }

    g_cnstr <- function(params) {
        design <- update(initial_design, params)
        user_cnstr <- evaluate(subject_to, design, optimization = TRUE)
        return(user_cnstr)
    }

    if (any(g_cnstr(tunable_parameters(initial_design)) > 0))
        warning("initial design is infeasible!")

    res <- nloptr::nloptr(
        x0 = tunable_parameters(initial_design),
        lb = tunable_parameters(lower_boundary_design),
        ub = tunable_parameters(upper_boundary_design),
        eval_f      = f_obj,
        eval_g_ineq = g_cnstr,
        opts        = opts,
        ...
    )

    if (res$status == 5 | res$status == 6)
        warning(res$message)

    return(list(
        design        = update(initial_design, res$solution),
        nloptr_return = res,
        call_args     = args
    ))

}
