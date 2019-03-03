#' Find optimal two-stage design by constraint minimization
#'
#' \code{minimize} takes an unconditioonal score
#' and a constraint set (or single constraint) of conditional and/or unconditional
#' scores and solves the corresponding constraint minimization problem
#' using \code{nloptr} (using COBYLA by default).
#'
#' @param objective objective function
#' @param subject_to constraint collection
#' @param initial_design initial guess (x0 for nloptr)
#' @param lower_boundary_design design specifying the lower boundary
#' @param upper_boundary_design design specifying the upper boundary
#' @param opts options list passed to nloptr
#' @param ... further optional arguments passed to \code{\link{nloptr}}
#'
#' @return \item{design}{ The resulting optimal design}
#'         \item{nloptr_return}{ Output of the corresponding nloptr call}
#'         \item{call_args}{ The arguments given to the optimization call}
#'
#'
#' @export
minimize <- function(
    objective,
    subject_to,
    initial_design,
    lower_boundary_design = NULL,
    upper_boundary_design = NULL,
    opts         =  list(
        algorithm   = "NLOPT_LN_COBYLA",
        xtol_rel    = 1e-5,
        maxeval     = 10000
    ),
    ...
) {

    args <- c(as.list(environment()), list(...))


    if(is(initial_design, "OneStageDesign")) {
        res <- do.call(.minimize_os, args = args)

    } else{
        if(is.null(lower_boundary_design) || is.null(upper_boundary_design)) {
            # compute boundaries if necessary

            cnstrs <- subject_to
            cnstrs@conditional_constraints <- list()

            bounds <- .starting_designs(objective, cnstrs, initial_design)

            if(is.null(lower_boundary_design))
                args$lower_boundary_design <- bounds$lb_design

            if(is.null(upper_boundary_design))
                args$upper_boundary_design <- bounds$ub_design

        }

        res <- do.call(.minimize, args = args)

    }

    return(res)
}



# minimization of two-stage designs
.minimize <- function(
    objective,
    subject_to,
    initial_design,
    lower_boundary_design = NULL,
    upper_boundary_design = NULL,
    opts,
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
                design@c1f - design@c1e + .1
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




# minimization of one-stage designs
.minimize_os <- function(
    objective,
    subject_to,
    initial_design,
    lower_boundary_design = NULL,
    upper_boundary_design = NULL,
    opts         =  list(
        algorithm   = "NLOPT_LN_COBYLA",
        xtol_rel    = 1e-5,
        maxeval     = 10000
    ),
    ...
) {

    if(is.null(lower_boundary_design))
        lower_boundary_design = OneStageDesign(5, 0)

    if(is.null(upper_boundary_design))
        upper_boundary_design = OneStageDesign(1000, 5)

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
            design@c1f - design@c1e
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




# compute boundary designs
.starting_designs <- function(
    objective,
    subject_to,
    initial_design,
    ...
) {
    os_design <- .minimize_os(objective, subject_to, OneStageDesign(200, 2))$design

    if(is(initial_design, "GroupSequentialDesign")) {
        lb_design <- GroupSequentialDesign(5,
                                           -1,
                                           os_design@c1f,
                                           2,
                                           rep(-2, length(initial_design@c2_pivots)))
        ub_design <- GroupSequentialDesign(os_design@n1,
                                           os_design@c1f,
                                           5,
                                           2*os_design@n1,
                                           rep(5, length(initial_design@c2_pivots)))
    } else{
        lb_design <- TwoStageDesign(5,
                                    -2,
                                    os_design@c1f,
                                    rep(2, length(initial_design@c2_pivots)),
                                    rep(-2, length(initial_design@c2_pivots)))
        ub_design <- TwoStageDesign(os_design@n1,
                                    os_design@c1f,
                                    5,
                                    rep(2*os_design@n1, length(initial_design@c2_pivots)),
                                    rep(5, length(initial_design@c2_pivots)))

    }

    return(list(
        lb_design = lb_design,
        ub_design = ub_design
    ))
}
