#' Two-stage designs
#'
#' \code{TwoStageDesign} is the fundamental design class of this package.
#' It consits of three first-stage parameters, two-stage two vectors
#' and additional parameters for implementation.
#'
#' \code{n1} denotes the sample size of the first stage.
#' \code{c1f} and \code{c1e} define decision boundary.
#' If the first-stage test statistic \code{Z_1} is larger than \code{c1e},
#' than the trial is stopped early for efficacy and the null hypothesis is
#' rejected.
#' If \code{Z_1} is smaller than \code{c1f}, than the null hypothesis is
#' accpected and the trial is stopped early for futility.
#' Otherwise, the trial enters in the second stage. In this case,
#' the stage-two sample size function \code{n2(Z_1)} and the stage-two
#' rejection boundary \code{c2(Z_1)} are computed.
#' As these are functions, they are approximated by a vector of pivots
#' denoted by \code{n2_pivots} and \code{c2_pivots}, respectively.
#' The slots \code{x1_norm_pivots} and \code{weights} are defined
#' for numerical implementation rules.
#' The generic implementation of this package is a Gaussian quadrature and
#' a corresponding design can be built by means of \link{gq_design}.
#' The user is free to implement own integration rules.
#' He has to be aware that all elements of \code{x1_norm_pivots} have
#' to be inside the interval [-1, 1] and have to be scaled by the
#' applied integration rule.
#' The parameter \code{tunable} is a logical one and indicates which
#' design parameters should be tunable for use in optimization.
#' Rounded sample size values can be applied by the logical parameter
#' \code{rounded}.
#'
#'
#' @slot n1 stage-one sample size
#' @slot c1f early stopping for futility boundary
#' @slot c1e early stopping for efficacy boundary
#' @slot n2_pivots vector of length order giving the values of n2 at the
#'     pivot points of the numeric integration rule [TODO: these are not available during construction]
#' @slot c2_pivots vector of length order giving the values of c2 at the
#'     pivot points of the numeric integration rule [TODO: these are not available during construction]
#' @slot x1_norm_pivots normalized pivots for integration rule (in [-1, 1])
#' @slot weights weights of conditional score values at x1_norm_pivots for
#'     approximating the integral over x1.
#' @slot tunable named logical vector indicating whether corresponding slot is considered a tunable parameter
#'
#' @exportClass TwoStageDesign
setClass("TwoStageDesign", representation(
        n1        = "numeric",
        c1f       = "numeric",
        c1e       = "numeric",
        n2_pivots = "numeric",
        c2_pivots = "numeric",
        x1_norm_pivots = "numeric",
        weights   = "numeric",
        tunable   = "logical"
    ))



#' @param ... optional arguments depending on implementation
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("TwoStageDesign", function(...) standardGeneric("TwoStageDesign"))


#' Create an object of class \code{TwoStageDesign}
#'
#' The argument \code{order} defines the order of the Gaussian quadrature rule.
#' If it is missing, the length of \code{n2_pivots} is used.
#'
#' If \code{order} is specified, \code{n2_pivots} and \code{c2_pivots} should
#' be one-dimensional. They are converted to vectors of length order
#' automatically.
#'
#' @param n1 cf. slot
#' @param c1f cf. slot
#' @param c1e cf. slot
#' @param n2_pivots cf. slot
#' @param c2_pivots cf. slot
#' @param order order (i.e. number of pivot points in the interior of [c1f, c1e])
#'     of the Gaussian quadrature rule to use for integration
#'
#' @return an object of class \code{TwoStageDesign}
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("TwoStageDesign", signature = "numeric",
function(n1, c1f, c1e, n2_pivots, c2_pivots, order = NULL, ...) {
    if(is.null(order)) {
        if(length(n2_pivots) != length(c2_pivots))
            stop("n2_pivots and c2_pivots must be of same length!")
        order <- length(n2_pivots)
    } else{
        n2_pivots <- rep(n2_pivots[1], order)
        c2_pivots <- rep(c2_pivots[1], order)
    }

    rule <- GaussLegendreRule(order)

    tunable <- logical(8) # initialize to all false
    tunable[1:5] <- TRUE
    names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")

    new("TwoStageDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
        c2_pivots = c2_pivots, x1_norm_pivots = rule$nodes, weights = rule$weights,
        tunable = tunable)
})




#' @rdname TwoStageDesign-class
#' @export
setGeneric("tunable_parameters", function(x, ...) standardGeneric("tunable_parameters"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("tunable_parameters", signature("TwoStageDesign"),
          function(x, ...) {
              res <- numeric(0)
              for (i in 1:length(x@tunable)) {
                  if (x@tunable[i])
                      res <- c(res, slot(x, names(x@tunable)[i]))
              }
              return(res)
          })


#' @rdname TwoStageDesign-class
#' @export
setGeneric("make_tunable", function(x, ...) standardGeneric("make_tunable"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("make_tunable", signature("TwoStageDesign"),
          function(x, ...) {
              params <- sapply(substitute(list(...))[-1], deparse)
              res    <- x
              for (i in 1:length(x@tunable)) {
                  if (names(x@tunable)[i] %in% params)
                      res@tunable[i] <- TRUE
              }
              return(res)
          })


#' @param x design
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("make_fixed", function(x, ...) standardGeneric("make_fixed"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("make_fixed", signature("TwoStageDesign"),
          function(x, ...) {
              params <- sapply(substitute(list(...))[-1], deparse)
              res    <- x
              for (i in 1:length(x@tunable)) {
                  if (names(x@tunable)[i] %in% params)
                      res@tunable[i] <- FALSE
              }
              return(res)
          })




#' @param params vector of design parameters (must be in same order as returned
#'     by \code{as.numeric(design)})
#' @param object object to update
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("update", signature("TwoStageDesign"),
    function(object, params, ...) {
        tunable_names <- names(object@tunable)[object@tunable]
        res <- object
        idx <- 1
        for (i in 1:length(tunable_names)) {
            slotname <- tunable_names[i]
            k <- length(slot(object, name = slotname))
            slot(res, name = slotname) <- params[idx:(idx + k - 1)]
            idx <- idx + k
        }
        return(res)
    })



#' @param d design object
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("n1", function(d, ...) standardGeneric("n1"))

#' @param round logical, should integer sample size or real sample size be
#'    returned?
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("n1", signature("TwoStageDesign"),
          function(d, round = TRUE, ...) {
              n1 <- d@n1
              if (round)
                  n1 <- round(n1)
              # if (round)
              #     n1 <- round(n1)
              return(n1)
          })



#' @param x1 stage-one outcome
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("n2", function(d, x1, ...) standardGeneric("n2"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("n2", signature("TwoStageDesign", "numeric"),
          function(d, x1, round = TRUE, ...) {
              res <- ifelse(x1 < d@c1f | x1 > d@c1e, 0, 1) *
                  pmax(
                      0,
                      stats::approx(
                          scaled_integration_pivots(d),
                          d@n2_pivots,
                          xout   = x1,
                          method = "linear",
                          rule   = 2
                      )$y
                  )
              if (round)
                  res <- round(res)
              return(res)
          })



#' @rdname TwoStageDesign-class
#' @export
setGeneric("n", function(d, x1, ...) standardGeneric("n"))

#' @describeIn TwoStageDesign overall sample size given stage-one outcome
#' @export
setMethod("n", signature("TwoStageDesign", "numeric"),
          function(d, x1, round = TRUE, ...) n2(d, x1, round, ...) + n1(d, round, ...))



#' @rdname TwoStageDesign-class
#' @export
setGeneric("c2", function(d, x1, ...) standardGeneric("c2"))

#' @rdname TwoStageDesign-class
#' @export
setMethod("c2", signature("TwoStageDesign", "numeric"),
          function(d, x1, ...) ifelse(x1 < d@c1f, Inf,
                                      ifelse(x1 > d@c1e, -Inf,
                                             stats::approx(scaled_integration_pivots(d), d@c2_pivots, xout = x1, method = "linear", rule = 2)$y)
          )
)


#' @rdname TwoStageDesign-class
#' @export
setGeneric("scaled_integration_pivots", function(d, ...) standardGeneric("scaled_integration_pivots"))

#' @describeIn TwoStageDesign get the actual pivots points (knots) of the numerical integration routine
#'     rule.
setMethod("scaled_integration_pivots", signature("TwoStageDesign"),
          function(d, ...){
              h <- (d@c1e - d@c1f) / 2
              return(h * d@x1_norm_pivots + (h + d@c1f))
          })




#' Show method for TwoStageDesign objects
#'
#' Only states the class itself. More information is given by the respective
#' \code{summary} method, cf. \link{TwoStageDesign-class}.
#'
#' @param object design to show
#'
#' @export
setMethod("show", signature(object = "TwoStageDesign"),
          function(object) cat("TwoStageDesign"))




#' Plot TwoStageDesign with optional set of conditional scores
#'
#' This function allows to plot the stage-two sample size and decision boundary
#' functions of a \code{TwoStageDesign} and
#' user-defined elements of the class \code{ConditionalScore}.
#'
#'
#' @param y not used
#' @param rounded should n-values be rounded?
#' @param k number of points to use for plotting
#'
#' @examples
#' order  <- 5L
#' design <- gq_design(50, 0, 2, rep(50.0, order), rep(2.0, order), order)
#' cp     <- ConditionalPower(dist = Normal(), prior = PointMassPrior(.4, 1))
#' plot(design, "Conditional Power" = cp)
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("plot", signature(x = "TwoStageDesign"),
          function(x, y = NULL, rounded = TRUE, ..., k = 100) {
              scores <- list(...)
              if (!all(sapply(scores, function(s) is(s, "ConditionalScore"))))
                  stop("optional arguments must be ConditionalScores")

              opts <- graphics::par(mfrow = c(1, length(scores) + 2))
              x1   <- seq(x@c1f, x@c1e, length.out = k)
              x2   <- seq(x@c1f - (x@c1e - x@c1f)/5, x@c1f - .01*(x@c1e - x@c1f)/5, length.out = k)
              x3   <- seq(x@c1e + .01*(x@c1e - x@c1f)/5, x@c1e + (x@c1e - x@c1f)/5, length.out = k)
              x4   <- seq(x@c1f - (x@c1e - x@c1f)/5, x@c1e + (x@c1e - x@c1f)/5, length.out = k)
              graphics::plot(x1, sapply(x1, function(z) n(x, z, round = rounded)), 'l',
                             xlim = c(min(x4), max(x4)),
                             ylim = c(0, 1.05 * max(sapply(x1, function(z) n(x, z, round = rounded)))),
                             main = "Overall sample size", ylab = "" , xlab = expression("x"[1]))
              graphics::lines(x2, sapply(x2, function(z) n(x, z, round = rounded)))
              graphics::lines(x3, sapply(x3, function(z) n(x, z, round = rounded)))
              graphics::plot(x4, c2(x, x4), 'l', main = "Stage-two critical value",
                             ylab = "", xlab = expression("x"[1]))
              if (length(scores) > 0) {
                  for (i in 1:length(scores)) {
                      graphics::plot(x1, evaluate(scores[[i]], x, x1), 'l',
                                     main = names(scores[i]), ylab = "", xlab = expression("x"[1]))
                      graphics::lines(x2, evaluate(scores[[i]], x, x2))
                      graphics::lines(x3, evaluate(scores[[i]], x, x3))
                  }
              }
              graphics::par(opts)
          })



#' Summarize TwoStageDesign objects with optional set of scores
#'
#' \code{summary} summarizes the first-stage of a \code{TwoStageDesign}
#' and objects of class \code{UnconditionalScore} that have to be
#' defined by the user.
#'
#' The first stage sample size and the two continuation decsion boundaries
#' are printed.
#' Furthermore, the user can define unconditional scores and these will be
#' evaluated and listed.
#'
#' @param object design object to summarize
#' @param rounded should rounded n-values be used?
#' @param ... optinal additional named UnconditionalScores
#'
#' @examples
#' order  <- 5L
#' design <- gq_design(50, 0, 2, rep(50.0, order), rep(2.0, order), order)
#' pow    <- integrate(ConditionalPower(dist = Normal(), prior = PointMassPrior(.4, 1)))
#' summary(design, "Power" = pow)
#'
#' @export
setMethod("summary", signature("TwoStageDesign"),
          function(object, ..., rounded = TRUE) {
              scores <- list(...)
              if (!all(sapply(scores, function(s) is(s, "UnconditionalScore"))))
                  stop("optional arguments must be UnconditionalScores")
              res <- list(
                  design = object,
                  scores = sapply(scores, function(s) evaluate(s, object, round = rounded, ...))
              )
              names(res$scores) <- names(scores)
              class(res) <- c("TwoStageDesignSummary", "list")
              return(res)
          })


#' Print object of class TwoStageDesignSummary
#'
#' @param x object to print
#' @param rounded should rounded n-values be used?
#' @param ... unused
#'
#' @export
print.TwoStageDesignSummary <- function(x, ..., rounded = TRUE) {
    cat("TwoStageDesign with:\n\r")
    cat(sprintf("     n1: %6.1f\n\r", n1(x$design, round = rounded)))
    cat(sprintf("    c1f: %6.1f\n\r", x$design@c1f))
    cat(sprintf("    c1e: %6.1f\n\r", x$design@c1e))
    if (length(x$scores) > 0) {
        cat("Unconditional scores:\n\r")
        for (i in 1:length(x$scores)) {
            cat(sprintf("    %s: %7.2f\n\r", names(x$scores)[i], x$scores[i]))
        }
        cat("\n\r")
    }
}



#' @describeIn TwoStageDesign simulate from the given design under parameter theta.
#'
#' @param nsim number of simulation runs
#' @param seed random seed
#' @param dist data distribution
#' @param theta location parameter of the data distribution
#'
#' @export
setMethod("simulate", signature("TwoStageDesign", "numeric"),
          function(object, nsim, dist, theta, seed = NULL, ...){
              if (!is.null(seed))
                  set.seed(seed)

              res <- data.frame(
                  theta  = rep(theta, nsim),
                  n1     = n1(object, round = TRUE),
                  c1f    = object@c1f,
                  c1e    = object@c1e
              )

              res$x1     <- simulate(dist, nsim = nsim, n = res$n1, theta = theta)
              res$n2     <- n2(object, res$x1, round = TRUE)
              res$c2     <- c2(object, res$x1)

              res$x2     <- simulate(dist, nsim = nsim, n = res$n2, theta = theta)
              res$reject <- res$x2 > res$c2 # check > vs. >=

              return(res)
          })
