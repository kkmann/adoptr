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
#' @param c2_monotone should the c2-function be forced to be monotoneously decreasing?
#' @param post_process should the sample sizes be integers?
#' @param opts options list passed to nloptr
#' @param ... further optional arguments passed to \code{\link{nloptr}}
#'
#' @export
minimize <- function(objective, subject_to, initial_design,
                     lower_boundary_design, upper_boundary_design,
                     c2_monotone = FALSE,
                     post_process = FALSE,
                     opts = list(
                         algorithm   = "NLOPT_LN_COBYLA",
                         xtol_rel    = 1e-5,
                         maxeval     = (200 * length(initial_design@c2_pivots)^2)
                     ), ...) {

        f_obj <- function(params) evaluate(objective, update(initial_design, params))

        g_cnstr <- function(params) {
            design <- update(initial_design, params)
            user_cnstr <- evaluate(subject_to, design)
            return(c(
                user_cnstr,
                design@c1f - design@c1e + ifelse( # ensure c1e > c1f if not one-stage
                    is(initial_design, "OneStageDesign"), 0, .1),
                if(c2_monotone == TRUE) diff(c2(design, scaled_integration_pivots(design))) # make c2() monotone if desired
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

        if(post_process == TRUE){
            n1 <- NULL
            n2_pivots <- NULL

            # Define continuous design as starting value and fix rounded sample sizes
            cont_design <- update(initial_design, res$solution)
            cont_design@n1 <- round(cont_design@n1)
            cont_design@n2_pivots <- round(cont_design@n2_pivots)
            cont_design <- make_fixed(cont_design, n1, n2_pivots)

            # Define new lower boundary design and fix rounded sample sizes
            lb_design <- lower_boundary_design
            lb_design@n1 <- cont_design@n1
            lb_design@n2_pivots <- cont_design@n2_pivots
            lb_design  <- make_fixed(lb_design, n1, n2_pivots)

            # Define new upper boundary design and fix rounded sample sizes
            ub_design <- upper_boundary_design
            ub_design@n1 <- cont_design@n1
            ub_design@n2_pivots <- cont_design@n2_pivots
            ub_design  <- make_fixed(ub_design, n1, n2_pivots)

            f_obj <- function(params) evaluate(objective, update(cont_design, params))

            g_cnstr <- function(params) {
                design <- update(cont_design, params)
                user_cnstr <- evaluate(subject_to, design)
                return(c(
                    user_cnstr,
                    design@c1f - design@c1e + ifelse( # ensure c1e > c1f if not one-stage
                        is(cont_design, "OneStageDesign"), 0, .1),
                    if(c2_monotone == TRUE) diff(c2(design, scaled_integration_pivots(design))) # make c2() monotone if desired
                ))
            }

            # Re-optimize c-values
            res2 <- nloptr::nloptr(
                x0 = tunable_parameters(cont_design),
                lb = tunable_parameters(lb_design),
                ub = tunable_parameters(ub_design),
                eval_f      = f_obj,
                eval_g_ineq = g_cnstr,
                opts = opts,
                ...
            )

            # Re-make parameters tunable for further use
            cont_design <- update(cont_design, res2$solution)
            cont_design <- make_tunable(cont_design, n1, n2_pivots)

            out <- list(
                "design" = cont_design,
                "nloptr_output" = res,
                "nloptr_output_post_processing" = res2
                )

            if(res$iterations == opts$maxeval){
                warning("Maximum number of iterations reached! Algorithm did probably not converge.")
            }

            return(out)

        } else{

            out <- list(
                "design" = update(initial_design, res$solution),
                "nloptr_output" = res
                )

            if(res$iterations == opts$maxeval){
                warning("Maximum number of iterations reached! Algorithm did probably not converge.")
            }


            return(out)
        }

        # TODO: error handling


    }
