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
#' @return a list with elements:
#'     \item{design}{ The resulting optimal design}
#'     \item{nloptr_return}{ Output of the corresponding nloptr call}
#'     \item{call_args}{ The arguments given to the optimization call}
#'
#' @examples
#' # Define Type one error rate
#' toer <- Power(Normal(), PointMassPrior(0.0, 1))
#'
#' # Define Power at delta = 0.4
#' pow <- Power(Normal(), PointMassPrior(0.4, 1))
#'
#' # Define expected sample size at delta = 0.4
#' ess <- ExpectedSampleSize(Normal(), PointMassPrior(0.4, 1))
#'
#' # Compute design minimizing ess subject to power and toer constraints
#' \dontrun{
#' minimize(
#'
#'    ess,
#'
#'    subject_to(
#'       toer <= 0.025,
#'       pow  >= 0.9
#'    ),
#'
#'    initial_design = TwoStageDesign(50, .0, 2.0, 60.0, 2.0, 5L)
#'
#' )
#' }
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

    res <- nloptr::nloptr(
        x0          = tunable_parameters(initial_design),
        lb          = tunable_parameters(lower_boundary_design),
        ub          = tunable_parameters(upper_boundary_design),
        eval_f      = f_obj,
        eval_g_ineq = g_cnstr,
        opts        = opts,
        ...
    )

    if (res$status == 5 | res$status == 6)
        warning(res$message)

    res <- list(
        design        = update(initial_design, res$solution),
        nloptr_return = res,
        call_args     = args
    )
    class(res) <- c("adoptrOptimizationResult", class(res))
    return(res)
}



#' @rawNamespace S3method(print, adoptrOptimizationResult)
print.adoptrOptimizationResult <- function(x, ...) {
    cat(design2str(x$design, TRUE))
}




#' Initial design
#'
#' The optimization method \code{\link{minimize}} requires an initial
#' design for optimization.
#' The function \code{get_initial_design} provides an initial guess based on a
#' fixed design that fulfills constraints on type I error rate and power.
#' Note that a situation-specific initial design may be much more efficient.
#'
#' @param theta the alternative effect size
#' @param alpha maximal type I error rate
#' @param beta maximale type II error rate
#' @param type is a two-stage, group-sequential, or one-stage design requried?
#' @param dist distribution of the test statistic
#' @param order desired integration order
#' @template dotdotdot
#'
#' @details
#' The distribution of the test statistic is specified by \code{dist}.
#' The default assumes a two-armed z-test.
#'
#' @examples
#' init <- get_initial_design(
#'    theta = 0.3,
#'    alpha = 0.025,
#'    beta  = 0.2,
#'    type  = "two-stage",
#'    dist  = Normal(two_armed = FALSE),
#'    order = 7L
#' )
#'
#' @export
get_initial_design <- function(theta, alpha, beta,
                               type = c("two-stage", "group-sequential", "one-stage"),
                               dist = Normal(), order = 7L, ...) {
    type <- match.arg(type)
    if (alpha <= 0 || alpha >= 1 || beta <= 0 || beta >= 1)
        stop("alpha and beta must be in (0, 1)!")
    theta <- ifelse(dist@two_armed, theta / sqrt(2), theta)
    if (is(dist, "Binomial")) {
        p_0   <- (theta + dist@rate_control + dist@rate_control) / 2
        theta <- theta / sqrt(p_0 * (1 - p_0))
    }
    c     <- quantile(dist, 1 - alpha, 1, 0)
    n     <- floor(2 * (c + quantile(dist, 1 - beta, 1, 0))^2 / theta^2)
    if (type == "one-stage")
        return(OneStageDesign(n, c))
    else if (type == "group-sequential")
        return(GroupSequentialDesign(n/2, 0, quantile(dist, 1 - alpha/2, 1, 0), n/2, c, order))
    else if (type == "two-stage")
        return(TwoStageDesign(n/2, 0, quantile(dist, 1 - alpha/2, 1, 0), n/2, c, order))
}





#' Boundary designs
#'
#' The optimization method \code{\link{minimize}} is based on the package
#' \code{nloptr}. This requires upper and lower boundaries for optimization.
#' Such boundaries can be computed via \code{lower_boundary_design}
#' respectively \code{upper_boundary_design}.
#' They are implemented by default in \code{\link{minimize}}.
#' Note that \code{\link{minimize}} allows the user to define its own
#' boundary designs, too.
#'
#' @param initial_design The initial design
#' @param n1 bound for the first-stage sample size n1
#' @param n2_pivots bound for the second-stage sample size n2
#' @param c1_buffer shift of the early-stopping boundaries from the initial ones
#' @param c2_buffer shift of the final decision boundary from the initial one
#' @param ... optional arguments
#'
#' The values \code{c1f} and \code{c1e} from the initial design are shifted
#' to \code{c1f - c1_buffer} and \code{c1e - c1_buffer} in
#' \code{get_lower_boundary_design}, respectively, to \cr
#' \code{c1f + c1_buffer} and \code{c1e + c1_buffer} in
#' \code{get_upper_boundary_design}.
#' This is handled analogously with \code{c2_pivots} and \code{c2_buffer}.
#'
#' @examples
#' initial_design <- TwoStageDesign(
#'   n1    = 25,
#'   c1f   = 0,
#'   c1e   = 2.5,
#'   n2    = 50,
#'   c2    = 1.96,
#'   order = 7L
#'   )
#' get_lower_boundary_design(initial_design)
#'
#' @rdname boundary-designs
#' @export
setGeneric("get_lower_boundary_design",
           function(initial_design, ...) standardGeneric("get_lower_boundary_design"))



#' @rdname boundary-designs
#' @export
setGeneric("get_upper_boundary_design",
           function(initial_design, ...) standardGeneric("get_upper_boundary_design"))



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("OneStageDesign"),
          function(initial_design, n1 = 1, c1_buffer = 2, ...) {
              lb_design <- OneStageDesign(n1, min(0, initial_design@c1f - c1_buffer))
              lb_design@tunable <- initial_design@tunable
              return(lb_design)

          })


#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("GroupSequentialDesign"),
          function(
              initial_design,
              n1        = 1,
              n2_pivots = 1,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              lb_design <- GroupSequentialDesign(
                  n1,
                  initial_design@c1f - c1_buffer,
                  initial_design@c1e - c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots - c2_buffer,
                  order = length(initial_design@c2_pivots)
              )
              lb_design@tunable <- initial_design@tunable
              return(lb_design)

          })



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("TwoStageDesign"),
          function(
              initial_design,
              n1        = 1,
              n2_pivots = 1,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              lb_design <- TwoStageDesign(
                  n1,
                  initial_design@c1f - c1_buffer,
                  initial_design@c1e - c1_buffer,
                  rep(n2_pivots, length(initial_design@c2_pivots)),
                  initial_design@c2_pivots - c2_buffer,
                  order = length(initial_design@c2_pivots)
              )
              lb_design@tunable <- initial_design@tunable
              return(lb_design)

          })





#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("OneStageDesign"),
          function(initial_design, n1 = 5 * initial_design@n1, c1_buffer = 2, ...) {
              ub_design <- OneStageDesign(n1, initial_design@c1f + c1_buffer)
              ub_design@tunable <- initial_design@tunable
              return(ub_design)
          })



#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("GroupSequentialDesign"),
          function(
              initial_design,
              n1        = 5 * initial_design@n1,
              n2_pivots = 5 * initial_design@n2_pivots,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              ub_design <- GroupSequentialDesign(
                  n1,
                  initial_design@c1f + c1_buffer,
                  initial_design@c1e + c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots + c2_buffer,
                  order = length(initial_design@c2_pivots)
              )
              ub_design@tunable <- initial_design@tunable
              return(ub_design)
          })


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("TwoStageDesign"),
          function(
              initial_design,
              n1        = 5 * initial_design@n1,
              n2_pivots = 5 * initial_design@n2_pivots,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              ub_design <- TwoStageDesign(
                  n1,
                  initial_design@c1f + c1_buffer,
                  initial_design@c1e + c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots + c2_buffer,
                  order = length(initial_design@c2_pivots)
              )
              ub_design@tunable <- initial_design@tunable
              return(ub_design)

          })
