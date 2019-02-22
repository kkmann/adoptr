.postprocess <- function(
    optimal_design,
    objective,
    subject_to,
    lower_boundary_design,
    upper_boundary_design,
    c2_monotone  = FALSE,
    opts         =  list(
        algorithm   = "NLOPT_LN_COBYLA",
        xtol_rel    = 1e-5,
        maxeval     = 10000
    ),
    ...
) {

    n1 <- NULL

    # if (is(initial_design, "OneStageDesign")) {
    #     # Define continuous design as starting value and fix rounded sample sizes
    #     cont_design <- update(initial_design, res$solution)
    #     cont_design@n1 <- round(cont_design@n1)
    #     cont_design <- make_fixed(cont_design, n1)
    #
    #     # Define new lower boundary design and fix rounded sample sizes
    #     lb_design <- update(cont_design, lower_boundary_design@c1f)
    #
    #     # Define new upper boundary design and fix rounded sample sizes
    #     ub_design <- update(cont_design, upper_boundary_design@c1f)
    #
    #
    # } else {
    # initial_design is not a one stage design
    n2_pivots <- NULL

    # define continuous design as starting value and fix rounded sample sizes
    post_design           <- optimal_design
    post_design@n1        <- n1(optimal_design, round = TRUE)
    post_design@n2_pivots <- round(optimal_design@n2_pivots)
    post_design           <- make_fixed(post_design, n1, n2_pivots)

    # define new lower boundary design and fix rounded sample sizes
    lb_design             <- lower_boundary_design
    lb_design@n1          <- post_design@n1
    lb_design@n2_pivots   <- post_design@n2_pivots
    lb_design             <- make_fixed(lb_design, n1, n2_pivots)

    # define new upper boundary design and fix rounded sample sizes
    ub_design             <- upper_boundary_design
    ub_design@n1          <- post_design@n1
    ub_design@n2_pivots   <- post_design@n2_pivots
    ub_design             <- make_fixed(ub_design, n1, n2_pivots)

    f_obj <- function(params) {
        evaluate(
            objective,
            update(post_design, params),
            optimization = TRUE
        )
    }

    g_cnstr <- function(params) {
        design <- update(post_design, params)
        cnstr  <- evaluate(subject_to, design, optimization = TRUE)
        return(c(
            cnstr,
            design@c1f - design@c1e + ifelse( # ensure c1e > c1f if not one-stage
                is(design, "OneStageDesign"), 0, .1),
            if (c2_monotone == TRUE) diff(c2(design, scaled_integration_pivots(design))) # make c2() monotone if desired
        ))
    }


    # re-optimize c-values
    res <- nloptr::nloptr(
        x0 = tunable_parameters(post_design),
        lb = tunable_parameters(lb_design),
        ub = tunable_parameters(ub_design),
        eval_f      = f_obj,
        eval_g_ineq = g_cnstr,
        opts        = opts,
        ...
    )

    if (res$status == 5 | res$status == 6)
        warning(res$message)

    # re-make parameters tunable for further use
    post_design <- update(post_design, res$solution)
    post_design <- make_tunable(post_design, n1, n2_pivots)

    # if (is(cont_design, "OneStageDesign")) {
    #     cont_design <- make_tunable(cont_design, n1)
    # } else{
    #     cont_design <- make_tunable(cont_design, n1, n2_pivots)
    # }

    return(list(
        design        = post_design,
        nloptr_return = res
    ))
}



#' Post-process an optimal design
#'
#' \code{postprocess} takes an optimal design and rounds its sample sizes.
#' The corresponding decision boundaries are re-computed such that the
#' constraints which were specified in the underlying \code{minimize()} call
#' are fulfilled.
#'
#' @param results An object obtained by \code{\link{minimize}}.
#'
#' @return \item{design}{ The resulting optimal design}
#'         \item{nloptr_return}{ Output of the corresponding nloptr call}
#'
#'
#' @export
postprocess <- function(results) {
    results$call_args$initial_design <- NULL
    do.call(.postprocess, args = c(list(optimal_design = results$design), results$call_args))
}
