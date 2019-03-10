#' Post-process decision boundaries of a two-stage design
#'
#' The results of a call to \code{\link{minimize}} is an design optimizing the
#' objective by relaxing all sample sizes to continuous variables.
#' In practice, sample sizes obviously need to be natural numbers, to
#' \code{postprocess} can be used to re-solve the same problem defined in
#' a call to minimize after rounding all sample sizes to integer values
#' (i.e., \code{n1} and \code{n2_pivots}).
#' These values are then keept constant during optimization and only the
#' decision boundaries are re-adjusted.
#' Note that this will only have substantial effect on small designs, where
#' the discretization error might be substantial.
#' Note tha even a non-preprocessed design will return rounded sample sizes by
#' default (cf. \code{\link{n}}).
#' The only difference is that the internal slots \code{n1} and \code{n2_pivots}
#' may still be real numbers and the corresponding critical values are not
#' adusted for the post-hoc discretization of these.
#'
#' @param results a list, expected to be the return value of a call to
#'   \code{\link{minimize}}.
#' @template dotdotdot
#'
#' @return list with two elements:
#'    \item{design}{the post-processed optimal design}
#'    \item{nloptr_return}{\code{\link[nloptr]{nloptr}} output}
#'
#' @seealso \code{\link{minimize}}
#'
#' @export
postprocess <- function(results, ...) {

    # calling arguments to minimize
    args <- results$call_args

    # change initial design to optimal design with integer sample sizes
    opt_design_int           <- results$design
    opt_design_int@n1        <- max(1, round(opt_design_int@n1))
    opt_design_int@n2_pivots <- pmax(0, round(opt_design_int@n2_pivots))

    # fix n1/n2 during optimization for design and boundaries
    n2_pivots      <- NULL # trick to fool R CMD check into accepting NSE below
    opt_design_int <- make_fixed(opt_design_int, n1, n2_pivots)
    lb_design      <- make_fixed(args$lower_boundary_design, n1, n2_pivots)
    ub_design      <- make_fixed(args$upper_boundary_design, n1, n2_pivots)

    # re-optimize with new constraints
    f_obj <- function(params) {
        evaluate(
            args$objective,
            update(opt_design_int, params),
            optimization = TRUE
        )
    }

    g_cnstr <- function(params) {
        design <- update(opt_design_int, params)
        return(evaluate(args$subject_to, design, optimization = TRUE))
    }

    res <- nloptr::nloptr(
        x0 = tunable_parameters(opt_design_int),
        lb = tunable_parameters(lb_design),
        ub = tunable_parameters(ub_design),
        eval_f      = f_obj,
        eval_g_ineq = g_cnstr,
        opts        = args$opts,
        ...
    )

    if (res$status == 5 | res$status == 6)
        warning(res$message)

    # reset al parameters to being tunable
    post_design <- update(opt_design_int, res$solution)
    post_design <- make_tunable(post_design, n1, n2_pivots)

    return(list(
        design        = post_design,
        nloptr_return = res
    ))


}
