#' Find optimal two-stage design by constraint minimization
#'
#' \code{minimize} takes an unconditional score and
#' a constraint set (or no constraint) and solves the corresponding
#' minimization problem using
#' \href{https://cran.r-project.org/package=nloptr}{\code{nloptr}}
#' (using COBYLA by default).
#' An initial design has to be defined. It is also possible to define
#' lower- and upper-boundary designs. If this is not done, the boundaries are
#' determined automatically heuristically.
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
#' @export
minimize <- function(
    objective,
    subject_to,
    initial_design,
    lower_boundary_design = get_lower_boundary_design(initial_design),
    upper_boundary_design = get_upper_boundary_design(initial_design),
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
        return(c(
            user_cnstr,
            design@c1f - design@c1e + ifelse( # ensure c1e > c1f if not one-stage
                is(initial_design, "OneStageDesign"), 0, .1)
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
