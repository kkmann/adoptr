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
#' @param check_constraints boolean whether constraints should be checked to be fulfilled
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
    check_constraints = FALSE,
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

    if(check_constraints){
        new_des <- update(initial_design, res$solution)
        for(constr in subject_to@unconditional_constraints){
            if(constr@rhs==0){
                if(evaluate(constr@score,new_des)-constr@rhs>=0.001){
                    warning(sprintf("The following constraint could not be fulfilled: %s (absolute tolerance: %s)", utils::capture.output(show(constr)), format(0.001)))
                }
            }
            else{
                if(evaluate(constr@score,new_des)-constr@rhs>=min(0.01*abs(constr@rhs),0.49)){
                    warning(sprintf("The following constraint could not be fulfilled: %s (relative tolerance: %s)", utils::capture.output(show(constr)), format(0.01)))
                }
            }

        }
        for(constr in subject_to@conditional_constraints){
            grid <- seq(new_des@c1f,new_des@c1e,length.out=10)
            if(constr@rhs==0){
                if(any(evaluate(constr@score,new_des,grid)-constr@rhs>=0.001)){
                    warning(sprintf("The following constraint could not be fulfilled: %s (absolute tolerance: %s)", utils::capture.output(show(constr)), format(0.001)))
                }
            }
            else{
                if(any(evaluate(constr@score,new_des,grid)-constr@rhs>=min(0.01*abs(constr@rhs),0.49))){
                    warning(sprintf("The following constraint could not be fulfilled: %s (relative tolerance: %s)", utils::capture.output(show(constr)), format(0.01)))
                }
            }

        }
    }

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
    cat(design2str(x$design, TRUE), "\n")
}




#' Initial design
#'
#' The optimization method \code{\link{minimize}} requires an initial
#' design for optimization.
#' This function provides a variety of possibilities to hand-craft designs that fulfill type I error and type II
#' error constraints which may be used as initial designs.
#'
#' @param theta the alternative effect size in the normal case, the
#' rate difference under the alternative in the binomial case
#' @param alpha maximal type I error rate
#' @param beta maximal type II error rate
#' @param type_design is a two-stage, group-sequential, or one-stage design required?
#' @param type_c2 either linear-decreasing c2-function according to inverse normal combination test or constant c2
#' @param type_n2 design of n2-function
#' @param dist distribution of the test statistic
#' @param cf first-stage futility boundary
#' @param ce first-stage efficacy boundary. Note that specifying this boundary implies that the type I error constraint might not be fulfilled anymore
#' @param info_ratio the ratio between first and second stage sample size
#' @param slope slope of n2 function
#' @param weight weight of first stage test statistics in inverse normal combination test
#' @param order desired integration order
#' @template dotdotdot
#'
#' @details
#' The distribution of the test statistic is specified by \code{dist}.
#' The default assumes a two-armed z-test.
#' The first stage efficacy boundary and the c2 boundary are chosen as Pocock-boundaries, so either ce=c2
#' if c2 is constant or ce=c, where the null hypothesis is rejected if w1 Z1+w2 Z2>c.
#' By specifying ce, it's clear that the boundaries are not Pocock-boundaries anymore, so the type I error
#' constraint may not be fulfilled.
#'
#' @examples
#' init <- get_initial_design(
#'    theta = 0.3,
#'    alpha = 0.025,
#'    beta  = 0.2,
#'    type_design="two-stage",
#'    type_c2="linear_decreasing",
#'    type_n2="linear_increasing",
#'    dist=Normal(),
#'    cf=0.7,
#'    info_ratio=0.5,
#'    slope=23,
#'    weight = 1/sqrt(3)
#' )
#'
#' @export
get_initial_design <- function(theta,
                                alpha,
                                beta,
                                type_design=c("two-stage","group-sequential","one-stage"),
                                type_c2=c("linear_decreasing","constant"),
                                type_n2=c("optimal","constant","linear_decreasing","linear_increasing"),
                                dist=Normal(),
                                cf=0,
                                ce,
                                info_ratio=0.5,
                                slope,
                                weight = sqrt(info_ratio),
                                order=7L,...){

    #match the inputs
    type_design <- match.arg(type_design)

    if(is(dist,"Student")){
        dist <- Normal(two_armed = dist@two_armed)
    }

    #raise warnings if c2 or n2 are mistakenly specified
    if(type_design=="group-sequential" || type_design == "one-stage"){
        if(type_design == "one-stage"){
            if(!missing(type_c2)){
                warning("The type of the c2 function is not relevant for one-stage designs.")
            }
        }
        if(!missing(type_n2)){
            warning("The type of the n2 function is not relevant for one-stage or group-sequential designs.")
        }}
    type_c2 <- match.arg(type_c2)
    type_n2 <- match.arg(type_n2)

    power <- 1-beta
    w1 <- weight
    w2 <- sqrt(1-w1**2)

    #check for valid inputs
    if (alpha <= 0 || alpha >= 1 || beta <= 0 || beta >= 1 || info_ratio <0 || info_ratio > 1 || weight >= 1 || weight <= 0){
        stop("alpha, beta, information ratio or weightsmust be in (0, 1)!")}

    #c2 and n2 type not relevant for one-stage designs, so we first finish that case
    if(type_design=="one-stage"){
        ce <- ifelse(is(dist,"Survival"),quantile(dist,1-alpha,1,1),quantile(dist,1-alpha,1,0))
        find_n <- function(n){
            (1-cumulative_distribution_function(dist,ce,n,theta))-power
        }
        n <- stats::uniroot(find_n, interval=c(0,1000),extendInt="upX")$root
        ifelse(is(dist,"Survival"),return(OneStageDesign(n,ce,event_rate=dist@event_rate)),
               return(OneStageDesign(n,ce)))
    }

    #case 1: constant c2 function
    if(type_c2=="constant"){
        #we assume ce=c2
        if(missing(ce)){
            find_c <- function(c){
                integrand_c <- function(z){
                    (1-stats::pnorm(c))*stats::dnorm(z)
                }
                (1-stats::pnorm(c))+stats::integrate(integrand_c,lower=cf,upper=c)$value-alpha
            }
            c2 <- stats::uniroot(find_c, interval=c(cf,5),extendInt="yes")$root
            ce <- c2
        }
        else{
            find_c2 <- function(c){
                integrand_c <- function(z){
                    (1-stats::pnorm(c))*stats::dnorm(z)
                }
                (1-stats::pnorm(ce))+stats::integrate(integrand_c,lower=cf,upper=ce)$value-alpha
            }
            c_try <- try(stats::uniroot(find_c2,interval=c(cf,5),extendInt="yes")$root,silent=TRUE)
            if("try-error" %in% class(c_try)){
                warning("Type I error constraint cannot be fulfilled.")
                find_c <- function(c){
                    integrand_c <- function(z){
                        (1-stats::pnorm(c))*stats::dnorm(z)
                    }
                    (1-stats::pnorm(c))+stats::integrate(integrand_c,lower=cf,upper=c)$value-alpha
                }
                c2 <- stats::uniroot(find_c, interval=c(cf,5),extendInt="yes")$root
            }
            else{
                c2 <- stats::uniroot(find_c2, interval=c(cf,5),extendInt="yes")$root
            }
        }
    }

    #case2: linear decreasing c2 function by inverse normal combination
    if(type_c2=="linear_decreasing"){
        #we assume ce=c*, where c* is the boundary for Z*=w_1*Z_1+w_2*Z_2
        if(missing(ce)){
            find_ce <- function(c){
                integrand_ce <- function(z){
                    (1-stats::pnorm((c-w1*z)/sqrt(1-w1**2)))*stats::dnorm(z)
                }
                (1-stats::pnorm(c))+stats::integrate(integrand_ce,lower=cf,upper=c)$value-alpha
            }
            c <- stats::uniroot(find_ce,interval=c(cf,5),extendInt="yes")$root
            ce <- c
        }
        else{
            find_ce <- function(c){
                integrand_ce <- function(z){
                    (1-stats::pnorm((c-w1*z)/sqrt(1-w1**2)))*stats::dnorm(z)
                }
                (1-stats::pnorm(ce))+stats::integrate(integrand_ce,lower=cf,upper=ce)$value-alpha
            }
            c_try <- try(stats::uniroot(find_ce,interval=c(cf,5),extendInt="yes")$root,silent=TRUE)
            if("try-error" %in% class(c_try)){
                warning("Type I error constraint cannot be fulfilled.")
                c <- ce
            }
            else{
                c <- stats::uniroot(find_ce,interval=c(cf,5),extendInt="yes")$root}
        }

        #we need to evaluate this function at each pivot
        critical_values <- function(z){
            (c-w1*z)/sqrt(1-w1**2)
        }
        oldnodes <- GaussLegendreRule(as.integer(order))$nodes
        h <- (ce - cf) / 2
        newnodes <- h * oldnodes + (h + cf)
        c2 <- sapply(newnodes, critical_values)


    }

    # for group-sequential designs, n2 is constant. We choose n so that the toer and tter are not exceeded
    if(type_design=="group-sequential"){
        if(type_c2=="constant"){
            find_n <- function(n){
                integrand_n <- function(z){
                    (1-cumulative_distribution_function(dist,c2,(1-info_ratio)*n,theta))*
                        probability_density_function(dist,z,(info_ratio)*n,theta)
                }
                (1-cumulative_distribution_function(dist,ce,info_ratio*n,theta))+stats::integrate(integrand_n,cf,ce)$value-power
            }
        }
        if(type_c2=="linear_decreasing"){
            find_n <- function(n){
                integrand_n <- function(z){
                    (1-cumulative_distribution_function(dist,(c-w1*z)/sqrt(1-w1**2),(1-info_ratio)*n,theta))*
                        probability_density_function(dist,z,(info_ratio)*n,theta)
                }
                (1-cumulative_distribution_function(dist,ce,info_ratio*n,theta))+stats::integrate(integrand_n,cf,ce)$value-power
            }
        }
        #choose n1 and n2 according to information ratio
        n <- stats ::uniroot(find_n,interval = c(0,1000),extendInt = "upX")$root
        n1 <- info_ratio*n
        n2 <- (1-info_ratio)*n
        #if(length(c2)!=1) n2 <- rep(n2,order)

        if(is(dist,"Survival")) design <- GroupSequentialDesign(n1,cf,ce,n2,c2,order=order, event_rate=dist@event_rate)
        else design <-  GroupSequentialDesign(n1,cf,ce,n2,c2,order=order)
        return(design)
    }

    #for two-stage designs, the four different n2-designs need to be considered
    if(type_design=="two-stage"){
        #constant n2 is the same as group-sequential, but returns a two-stage instead of a group sequential design
        if(type_n2=="constant"){
            if(type_c2=="constant"){
                find_n <- function(n){
                    integrand_n <- function(z){
                        (1-cumulative_distribution_function(dist,c2,(1-info_ratio)*n,theta))*
                            probability_density_function(dist,z,(info_ratio)*n,theta)
                    }
                    (1-cumulative_distribution_function(dist,ce,info_ratio*n,theta))+stats::integrate(integrand_n,cf,ce)$value-power
                }
            }
            if(type_c2=="linear_decreasing"){
                find_n <- function(n){
                    integrand_n <- function(z){
                        (1-cumulative_distribution_function(dist,(c-w1*z)/sqrt(1-w1**2),(1-info_ratio)*n,theta))*
                            probability_density_function(dist,z,(info_ratio)*n,theta)
                    }
                    (1-cumulative_distribution_function(dist,ce,info_ratio*n,theta))+stats::integrate(integrand_n,cf,ce)$value-power
                }
            }
            #choose n1 and n2 according to information ratio
            n <- stats ::uniroot(find_n,interval = c(0,1000),extendInt = "upX")$root
            n1 <- info_ratio*n
            n2 <- (1-info_ratio)*n
            if(length(c2)!=1) n2 <- rep(n2,order)


            if(is(dist,"Survival")) design <- TwoStageDesign(n1,cf,ce,n2,c2,order=order, event_rate=dist@event_rate)
            else design <-  TwoStageDesign(n1,cf,ce,n2,c2,order=order)
            return(design)
        }

        #the optimal n2 function is evaluated at each pivot and therefore, it is not constant. The first stage
        #sample size is determined like before.
        if(type_n2=="optimal" || type_n2=="linear_increasing" || type_n2=="linear_decreasing"){
            if(length(c2)==1) c2 <- rep(c2,order)

            if(type_c2=="constant"){
                find_n <- function(n){
                    integrand_n <- function(z){
                        (1-cumulative_distribution_function(dist,c2,(1-info_ratio)*n,theta))*
                            probability_density_function(dist,z,(info_ratio)*n,theta)
                    }
                    (1-cumulative_distribution_function(dist,ce,info_ratio*n,theta))+stats::integrate(integrand_n,cf,ce)$value-power
                }
            }

            if(type_c2=="linear_decreasing"){
                #first, find the n1 sample size. We use the same sample size as we have under constant n2-function
                find_n <- function(n){
                    integrand_n <- function(z){
                        (1-cumulative_distribution_function(dist,(c-w1*z)/sqrt(1-w1**2),(1-info_ratio)*n,theta))*
                            probability_density_function(dist,z,(info_ratio)*n,theta)
                    }
                    (1-cumulative_distribution_function(dist,ce,info_ratio*n,theta))+stats::integrate(integrand_n,cf,ce)$value-power
                }
            }
            n <- stats ::uniroot(find_n,interval = c(0,1000),extendInt = "upX")$root
            n1 <- info_ratio*n

            #find the right value of the n2 function at each pivot
            find_n2 <- function(n2,pivot){
                integrand_n2 <- function(z){
                    (1-cumulative_distribution_function(dist,pivot,n2,theta))*
                        probability_density_function(dist,z,n1,theta)
                }
                (1-cumulative_distribution_function(dist,ce,n1,theta))+
                    stats::integrate(integrand_n2,cf,ce)$value-power
            }
            n2 <- rep(0.0,order)
            for(i in (1:order)){
                n_try <- suppressWarnings(try(stats::uniroot(f=find_n2,interval=c(0,1000),extendInt="upX",
                                                             pivot=c2[i])$root,silent=TRUE))
                if("try-error" %in% class(n_try) & type_n2 == "optimal"){
                    stop("Optimal n2 function cannot be calculated. Please reduce efficacy boundary or the information ratio.") # nocov
                }
                if("try-error" %in% class(n_try) & (type_n2 == "linear_decreasing" || type_n2 == "linear_increasing") & missing(slope)){
                    stop("Please specify a slope or reduce efficacy boundary or the information ratio.") # nocov
                }
                if("try-error" %in% class(n_try) & (type_n2 == "linear_decreasing" || type_n2 == "linear_increasing") & !missing(slope)){
                    n2[i] <- 0.1 # nocov
                }
                else{
                    n2[i] <- stats::uniroot(f=find_n2,interval=c(0,1000),extendInt="upX",
                                            pivot=c2[i])$root
                }
            }

            if(type_n2=="optimal"){
                if(is(dist,"Survival")) design <- TwoStageDesign(n1,cf,ce,n2,c2, event_rate=dist@event_rate)
                else design <-  TwoStageDesign(n1,cf,ce,n2,c2)
                return(design)
            }

            if(type_n2=="linear_decreasing"){
                interim_design <- TwoStageDesign(n1,cf,ce,n2,c2)
                if(missing(slope)){
                    #find the minimum and maximal values of n2
                    n2_min <- n2(interim_design,ce)
                    n2_max <- n2(interim_design,cf)

                    #compute the slope of n2 function
                    slope <- (n2_min-n2_max)/(ce-cf)
                    if(slope==0){
                        slope <- -n2_max/(ce-cf) # nocov
                    }
                }

                if(slope>0){
                    stop("For linear_decreasing, slope needs to be smaller zero")
                }
                start_value <- -slope*ce+.1

                #find y-intercept of n2 function
                find_y_intercept <- function(b,slope){
                    integrand_b <- function(z){
                        (1-cumulative_distribution_function(dist,c2(interim_design,z),slope*z+b,theta))*
                            probability_density_function(dist,z,n1,theta)
                    }
                    (1-cumulative_distribution_function(dist,ce,n1,theta))+
                        stats::integrate(integrand_b,cf,ce)$value-power
                }

                t <- suppressWarnings(try(stats::uniroot(find_y_intercept,interval=c(start_value,1000),extendInt = "upX",slope=slope)$root,silent=TRUE))
                if("try-error" %in% class(t)){
                    warning("Absolute value of slope too high. It was automatically reduced.")
                }
                while(("try-error" %in% class(t))&&(slope<=0)){
                    slope <- slope+1
                    start_value <- -slope*ce+.1
                    t <- suppressWarnings(try(stats::uniroot(find_y_intercept,interval=c(start_value,1000),extendInt = "upX",slope=slope)$root,silent=TRUE))
                }
                if(slope>0){# nocov start
                    if(is(dist,"Survival")) design <- TwoStageDesign(n1,cf,ce,n2,c2, event_rate=dist@event_rate)
                    else design <-  TwoStageDesign(n1,cf,ce,n2,c2)
                    return(design)# nocov end
                }

                y_intercept <- stats::uniroot(find_y_intercept,interval=c(start_value,1000),extendInt = "upX",slope=slope)$root

                n2func <- function(z){
                    slope*z+y_intercept
                }

                oldnodes <- GaussLegendreRule(as.integer(order))$nodes
                h <- (ce - cf) / 2
                newnodes <- h * oldnodes + (h + cf)

                #evaluate n2 function at pivots
                n2 <- n2func(newnodes)
                if(is(dist,"Survival")) design <- TwoStageDesign(n1,cf,ce,n2,c2, event_rate=dist@event_rate)
                else design <-  TwoStageDesign(n1,cf,ce,n2,c2)
                return(design)
            }

            if(type_n2=="linear_increasing"){
                #find the minimum and maximal values of n2
                interim_design <- TwoStageDesign(n1,cf,ce,n2,c2)
                if(missing(slope)){
                    #find the minimum and maximal values of n2
                    n2_min <- n2(interim_design,ce)
                    n2_max <- n2(interim_design,cf)


                    #compute the slope of n2 function
                    slope <- (n2_max-n2_min)/(ce-cf)
                    if(slope==0){
                        slope <- n2_max/(ce-cf) # nocov
                    }
                }
                if(slope<0){
                    stop("For linear_increasing, slope needs to be greater zero")
                }


                #find y-intercept of n2 function
                find_y_intercept <- function(b,slope){
                    integrand_b <- function(z){
                        (1-cumulative_distribution_function(dist,c2(interim_design,z),slope*z+b,theta))*
                            probability_density_function(dist,z,n1,theta)
                    }
                    (1-cumulative_distribution_function(dist,ce,n1,theta))+
                        stats::integrate(integrand_b,cf,ce)$value-power
                }

                start_value <- -slope*cf+.1

                #it is possible that the n2 function gets negative what leads to errors,
                #so we need to reduce the slope to be positive
                t <- suppressWarnings(try(stats::uniroot(find_y_intercept,interval=c(start_value,1000),extendInt = "upX",slope=slope)$root,silent=TRUE))

                if("try-error" %in% class(t)){
                    warning("Absolute value of slope too high. It was automatically reduced.")
                }
                while(("try-error" %in% class(t))&&(slope>=0)){
                    slope <- slope-1
                    start_value <- -slope*cf+.1
                    t <- suppressWarnings(try(stats::uniroot(find_y_intercept,interval=c(start_value,1000),extendInt = "upX",slope=slope)$root,silent=TRUE))
                }

                if(slope<0){# nocov start
                    if(is(dist,"Survival")) design <- TwoStageDesign(n1,cf,ce,n2,c2, event_rate=dist@event_rate)
                    else design <-  TwoStageDesign(n1,cf,ce,n2,c2)
                    return(design)# nocov end
                }

                y_intercept <- stats::uniroot(find_y_intercept,interval=c(start_value,1000),extendInt = "upX",slope=slope)$root

                n2func <- function(z){
                    slope*z+y_intercept
                }

                oldnodes <- GaussLegendreRule(as.integer(order))$nodes
                h <- (ce - cf) / 2
                newnodes <- h * oldnodes + (h + cf)

                #evaluate n2 function at pivots
                n2 <- n2func(newnodes)
                if(is(dist,"Survival")) design <- TwoStageDesign(n1,cf,ce,n2,c2, event_rate=dist@event_rate)
                else design <-  TwoStageDesign(n1,cf,ce,n2,c2)
                return(design)
            }

        }
    }
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
