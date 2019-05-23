#' Two-stage designs
#'
#' \code{TwoStageDesign} is the fundamental design class of the
#' \pkg{\link{adoptr}} package.
#' Formally, we represent a generic two-stage design as a five-tuple
#' \ifelse{html}{\out{(n<sub>1</sub>, c<sub>1</sub><sup>f</sup>, c<sub>1</sub><sup>e</sup>, n<sub>2</sub>(&middot;), c<sub>2</sub>(&middot;))}}{\eqn{\big(n_1, c_1^f, c_1^e, n_2(\cdot), c_2(\cdot)\big)}}.
#' Here, \ifelse{html}{\out{n<sub>1</sub>}}{\eqn{n_1}} is the first-stage sample
#' size (per group), \ifelse{html}{\out{c<sub>1</sub><sup>f</sup>}}{\eqn{c_1^f}}
#' and \ifelse{html}{\out{c<sub>1</sub><sup>e</sup>}}{\eqn{c_1^e}} are
#' boundaries for early stopping for futility and efficacy, respectively.
#' Since the trial design is a two-stage design, the elements
#' \ifelse{html}{\out{n<sub>2</sub>(&middot;)}}{\eqn{n_2(\cdot)}} (stage-two sample
#' size) and \ifelse{html}{\out{c<sub>2</sub>(&middot;)}}{\eqn{c_2(\cdot)}}
#' (stage-two critical value) are functions of the first-stage outcome
#' \ifelse{html}{\out{X<sub>1</sub>=x<sub>1</sub>}}{\eqn{X_1=x_1}}.
#' \ifelse{html}{\out{X<sub>1</sub>}}{\eqn{X_1}} denotes the first-stage test
#' statistic. A brief description on this definition of two-stage designs can be
#' read \href{https://kkmann.github.io/adoptr/articles/adoptr.html}{here}.
#' For available methods, see the 'See Also' section at the end of this page.
#'
#' @slot n1 cf. parameter 'n1'
#' @slot c1f cf. parameter 'c1f'
#' @slot c1e cf. parameter 'c1e'
#' @slot n2_pivots vector of length 'order' giving the values of n2 at the
#'     pivot points of the numeric integration rule
#' @slot c2_pivots vector of length order giving the values of c2 at the
#'     pivot points of the numeric integration rule
#' @slot x1_norm_pivots normalized pivots for integration rule (in [-1, 1])
#'     the actual pivots are scaled to the interval [c1f, c1e] and can be
#'     obtained by the internal method \cr
#'     \code{adoptr:::scaled_integration_pivots(design)}
#' @slot weights weights of of integration rule at \code{x1_norm_pivots} for
#'     approximating integrals over \code{x1}
#' @slot tunable named logical vector indicating whether corresponding slot is
#'     considered a tunable parameter (i.e. whether it can be changed during
#'     optimization via \code{\link{minimize}} or not; cf. \cr
#'     \code{\link{make_fixed}})
#'
#' @seealso For accessing sample sizes and critical values safely, see methods in
#' \code{\link{n}} and \code{\link{c2}}; for modifying behaviour during optimizaton
#' see \code{\link{make_tunable}}; to convert between S4 class represenation and
#' numeric vector, see \code{\link{tunable_parameters}}; for simulating from a given
#' design, see \code{\link[=simulate,TwoStageDesign,numeric-method]{simulate}};
#' for plotting see \code{\link{plot,TwoStageDesign-method}}.
#' Both \link[=GroupSequentialDesign-class]{group-sequential} and
#' \link[=OneStageDesign]{one-stage designs} (!) are implemented as subclasses of
#' \code{TwoStageDesign}.
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

#' @param n1 stage-one sample size
#'
#' @rdname TwoStageDesign-class
#' @export
setGeneric("TwoStageDesign", function(n1, ...) standardGeneric("TwoStageDesign"))

#' @template c1f
#' @template c1e
#' @param n2_pivots numeric vector, stage-two sample size on the integration
#' pivot points
#' @param c2_pivots numeric vector, stage-two critical values on the integration
#' pivot points
#' @template order
#' @template dotdotdot
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("TwoStageDesign", signature = "numeric",
    function(n1, c1f, c1e, n2_pivots, c2_pivots, order = NULL, ...) {

        if (length(n2_pivots) != length(c2_pivots))
            stop("n2_pivots and c2_pivots must be of same length!")
        if (is.null(order)) {
            order <- length(n2_pivots)
        } else if (length(n2_pivots) != order) {
            n2_pivots <- rep(n2_pivots[1], order)
            c2_pivots <- rep(c2_pivots[1], order)
        }

        rule           <- GaussLegendreRule(as.integer(order))
        tunable        <- logical(8) # initialize to all false
        tunable[1:5]   <- TRUE
        names(tunable) <- c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots", "x1_norm_pivots", "weights", "tunable")

        new("TwoStageDesign", n1 = n1, c1f = c1f, c1e = c1e, n2_pivots = n2_pivots,
            c2_pivots = c2_pivots, x1_norm_pivots = rule$nodes, weights = rule$weights,
            tunable = tunable)

    })




#' Switch between numeric and S4 class representation of a design
#'
#' Get tunable parameters of a design as numeric vector via
#' \code{tunable_parameters} or \code{update} a design object with a suitable
#' vector of values for its tunable parameters.
#'
#' @param    object    \code{TwoStageDesign} object to update
#' @template dotdotdot
#'
#' @details
#' The \code{tunable} slot of a \code{\link{TwoStageDesign}} stores information about
#' the set of design parameters which are considered fixed (not changed during
#' optimization) or tunable (changed during optimization).
#' For details on how to fix certain parameters or how to make them tunable
#' again, see \code{\link{make_fixed}} and \code{\link{make_tunable}}.
#'
#' @examples
#' design  <- TwoStageDesign(25, 0, 2, 25, 2, order = 5)
#' tunable_parameters(design)
#' design2 <-update(design, tunable_parameters(design) + 1)
#' tunable_parameters(design2)
#'
#' @seealso \code{\link{TwoStageDesign}}
#'
#' @export
setGeneric("tunable_parameters", function(object, ...) standardGeneric("tunable_parameters"))

#' @rdname tunable_parameters
#' @export
setMethod("tunable_parameters", signature("TwoStageDesign"),
          function(object, ...) {
              res <- numeric(0)
              for (i in 1:length(object@tunable)) {
                  if (object@tunable[i])
                      res <- c(res, slot(object, names(object@tunable)[i]))
              }
              return(res)
          })

#' @param params vector of design parameters, must be in same order as returned
#'   by \cr
#'   \code{tunable_parameters}
#'
#' @rdname tunable_parameters
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



#' Fix parameters during optimization
#'
#' The methods \code{make_fixed} and \code{make_tunable} can be used to modify
#' the 'tunability' status of parameters in a \code{\link{TwoStageDesign}}
#' object.
#' Tunable parameters are optimized over, non-tunable ('fixed') parameters are
#' considered given and not altered during optimization.
#'
#' @param x \code{TwoStageDesign} object
#' @param ... unquoted names of slots for which the tunability status should be
#' changed.
#'
#' @examples
#' design <- TwoStageDesign(25, 0, 2, 25, 2, order = 5)
#' # default: all parameters are tunable (except integration pivots,
#' # weights and tunability status itself)
#' design@tunable
#'
#' # make n1 and the pivots of n2 fixed (not changed during optimization)
#' design <- make_fixed(design, n1, n2_pivots)
#' design@tunable
#'
#' # make them tunable again
#' design <- make_tunable(design, n1, n2_pivots)
#' design@tunable
#'
#' @seealso \code{\link{TwoStageDesign}}, \code{\link{tunable_parameters}} for
#' converting tunable parameters of a design object to a numeric vector (and back),
#' and \code{\link{minimize}} for the actual minimzation procedure
#'
#' @export
setGeneric("make_tunable", function(x, ...) standardGeneric("make_tunable"))

#' @rdname make_tunable
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



#' @rdname make_tunable
#' @export
setGeneric("make_fixed", function(x, ...) standardGeneric("make_fixed"))

#' @rdname make_tunable
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









#' @rdname n
#' @export
setGeneric("n1", function(d, ...) standardGeneric("n1"))

#' @rdname n
#' @export
setMethod("n1", signature("TwoStageDesign"),
          function(d, round = TRUE, ...) {
              n1 <- d@n1
              if (round)
                  n1 <- round(n1)
              return(n1)
          })



#' @rdname n
#' @export
setGeneric("n2", function(d, x1, ...) standardGeneric("n2"))

#' @rdname n
#' @export
setMethod("n2", signature("TwoStageDesign", "numeric"),
          function(d, x1, round = TRUE, ...) {
              res <- ifelse(x1 < d@c1f | x1 > d@c1e, 0, 1) *
                  pmax(
                      0,
                      stats::splinefun(
                          scaled_integration_pivots(d),
                          d@n2_pivots,
                          method = "monoH.FC"
                        )(x1)
                  )
              if (round)
                  res <- round(res)
              return(res)
          })



#' Query sample size of a design
#'
#' Methods to access the stage-one, stage-two, or overall sample size of a
#' \code{\link{TwoStageDesign}}.
#' \code{n1} returns the first-stage sample size of a design,
#' \code{n2} the stage-two sample size conditional on the stage-one test
#' statistic and \code{n} the overall sample size \code{n1 + n2}.
#' Internally, objects of the class \code{TwoStageDesign} allow non-natural,
#' real sample sizes to allow smooth optimization (cf. \code{\link{minimize}} for
#' details).
#' The optional argument \code{round} allows to switch between the internal
#' real representation and a rounded version (rounding to the next positive
#' integer).
#'
#' @examples
#' design <- TwoStageDesign(
#'    n1    = 25,
#'    c1f   = 0,
#'    c1e   = 2.5,
#'    n2    = 50,
#'    c2    = 1.96,
#'    order = 7L
#' )
#'
#' n1(design) # 25
#' design@n1 # 25
#'
#' n(design, x1 = 2.2) # 75
#'
#'
#' @template d
#' @template x1
#' @template round
#' @template dotdotdot
#'
#' @seealso \code{\link{TwoStageDesign}}, see \code{\link{c2}} for accessing
#' the critical values
#'
#' @rdname n
#' @export
setGeneric("n", function(d, x1, ...) standardGeneric("n"))

#' @rdname n
#' @export
setMethod("n", signature("TwoStageDesign", "numeric"),
          function(d, x1, round = TRUE, ...) n2(d, x1, round, ...) + n1(d, round, ...))





#' Query critical values of a design
#'
#' Methods to access the stage-two critical values of a
#' \code{\link{TwoStageDesign}}.
#' \code{c2} returns the stage-two critical value conditional on the stage-one test
#' statistic.
#'
#' @examples
#' design <- TwoStageDesign(
#'   n1    = 25,
#'   c1f   = 0,
#'   c1e   = 2.5,
#'   n2    = 50,
#'   c2    = 1.96,
#'   order = 7L
#' )
#'
#' c2(design, 2.2) # 1.96
#' c2(design, 3.0) # -Inf
#' c2(design, -1.0) # Inf
#'
#' @template d
#' @template x1
#' @template dotdotdot
#'
#' @seealso \code{\link{TwoStageDesign}}, see \code{\link{n}} for accessing
#' the sample size of a design
#'
#' @examples
#' design <- TwoStageDesign(
#'    n1    = 25,
#'    c1f   = 0,
#'    c1e   = 2.5,
#'    n2    = 50,
#'    c2    = 1.96,
#'    order = 7L
#' )
#'
#' c2(design, 2.2) # 1.96
#' c2(design, 3.0) # -Inf
#' c2(design, -1.0) # Inf
#'
#' @rdname critical-values
#' @export
setGeneric("c2", function(d, x1, ...) standardGeneric("c2"))

#' @rdname critical-values
#' @export
setMethod("c2", signature("TwoStageDesign", "numeric"),
          function(d, x1, ...) ifelse(x1 < d@c1f, Inf,
                                      ifelse(x1 > d@c1e, -Inf,
                                             stats::splinefun(
                                                 scaled_integration_pivots(d),
                                                 d@c2_pivots,
                                                 method = "monoH.FC"
                                                 )(x1)
                                             )
          )
)





# internal, get integration pivots scales to [c1f, c1e]
setGeneric("scaled_integration_pivots", function(d, ...) standardGeneric("scaled_integration_pivots"))

setMethod("scaled_integration_pivots", signature("TwoStageDesign"),
          function(d, ...){
              h <- (d@c1e - d@c1f) / 2
              return(h * d@x1_norm_pivots + (h + d@c1f))
          })





#' @param object design to show or summarize
#'
#' @rdname TwoStageDesign-class
#' @export
setMethod("show", signature(object = "TwoStageDesign"),
          function(object) cat(class(object)[1]))




#' Plot \code{TwoStageDesign} with optional set of conditional scores
#'
#' This method allows to plot the stage-two sample size and decision boundary
#' functions of a chosen design.
#'
#' \code{\link{TwoStageDesign}} and
#' user-defined elements of the class \code{\link[=Scores]{ConditionalScore}}.
#'
#' @template plot
#' @param rounded should n-values be rounded?
#' @param k number of points to use for plotting
#' @param ... further named \code{ConditinonalScores} to plot for the design
#'
#' @seealso \code{\link{TwoStageDesign}}
#'
#' @examples
#' design <- TwoStageDesign(50, 0, 2, 50, 2, 5)
#' cp     <- ConditionalPower(dist = Normal(), prior = PointMassPrior(.4, 1))
#' plot(design, "Conditional Power" = cp)
#'
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
                      y <- list(
                          left =   evaluate(scores[[i]], x, x2),
                          middle = evaluate(scores[[i]], x, x1),
                          right =  evaluate(scores[[i]], x, x3)
                      )
                      expand <- .05*(max(do.call(c, y)) - min(do.call(c, y)))
                      graphics::plot(
                          x1, y$middle,
                          'l',
                          xlim = c(min(x4), max(x4)),
                          ylim = c(min(do.call(c, y)) - expand, max(do.call(c, y)) + expand),
                          main = names(scores[i]),
                          ylab = "",
                          xlab = expression("x"[1])
                      )
                      graphics::lines(x2, y$left)
                      graphics::lines(x3, y$right)
                  }
              }
              graphics::par(opts)
          })



#' @details
#' \code{summary} can be used to quickly compute and display basic facts about
#' a TwoStageDesign.
#' An arbitrary number of names \code{\link[=Scores]{UnconditionalScore}} objects can be
#' provided via the optional arguments \code{...} and are included in the summary displayed using
#' \code{\link{print}}.
#'
#' @param rounded should rounded n-values be used?
#'
#' @examples
#' design <- TwoStageDesign(50, 0, 2, 50.0, 2.0, 5)
#' pow    <- Power(Normal(), PointMassPrior(.4, 1))
#' summary(design, "Power" = pow)
#'
#' @rdname TwoStageDesign-class
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

#' @param  x return value of call to \code{summary}
#' @template round
#'
#' @rdname TwoStageDesign-class
#' @export
print.TwoStageDesignSummary <- function(x, ..., round = TRUE) {
    x1 <- seq(x$design@c1f, x$design@c1e, length.out = 1000)
            cat(sprintf("%s with:\n\r", class(x$design)[1]))
    cat(sprintf("     n1: %6.2f\n\r", n1(x$design, round)))
    cat(sprintf("    c1f: %6.2f\n\r", x$design@c1f))
    cat(sprintf("    c1e: %6.2f\n\r", x$design@c1e))
    cat(sprintf(" max n2: %6.2f\n\r", max(n2(x$design, x1, round))))
    cat(sprintf(" min n2: %6.2f\n\r", min(n2(x$design, x1, round))))
    cat(sprintf("%i integration pivots at: ", length(x$design@x1_norm_pivots)))
    cat(paste0(sprintf("%.2f", scaled_integration_pivots(x$design)), collapse = ", "))
    cat("\n\r    integration weights: ")
    cat(paste0(sprintf("%.2f", x$design@weights), collapse = ", "))
    cat("\n\r")
    if (length(x$scores) > 0) {
        cat("Unconditional scores:\n\r")
        for (i in 1:length(x$scores)) {
            cat(sprintf("    %10s: %7.3f\n\r", names(x$scores)[i], x$scores[i]))
        }
        cat("\n\r")
    }
}



#' Draw samples from a two-stage design
#'
#' \code{simulate} allows to draw samples from a given
#' \code{\link{TwoStageDesign}}.
#'
#' @param object \code{TwoStageDesign} to draw samples from
#' @param nsim number of simulation runs
#' @param seed random seed
#' @param dist data distribution
#' @param theta location parameter of the data distribution
#' @template dotdotdot
#'
#' @return \code{simulate()} returns a \code{data.frame} with \code{nsim}
#' rows and for each row (each simulation run) the following columns \itemize{
#' \item{theta }{The effect size}
#' \item{n1 }{First-stage sample size}
#' \item{c1f }{Stopping for futility boundary}
#' \item{c1e }{Stopping for efficacy boundary}
#' \item{x1 }{First-stage outcome}
#' \item{n2 }{Resulting second-stage sample size after observing x1}
#' \item{c2 }{Resulting second-stage decision-boundary after observing x1}
#' \item{x2 }{Second-stage outcome}
#' \item{reject }{Decision whether the null hypothesis is rejected or not}
#' }
#'
#' @examples
#' design <- TwoStageDesign(25, 0, 2, 25, 2, order = 5)
#' # draw samples assuming two-armed design
#' simulate(design, 10, Normal(), .3, 42)
#'
#' @seealso \code{\link{TwoStageDesign}}
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
