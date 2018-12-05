#' @rdname Design-class
#' @export
setGeneric("early_stopping_bounds", function(d, ...) standardGeneric("early_stopping_bounds"))


#' @rdname Design-class
#' @export
setGeneric("n1", function(d, ...) standardGeneric("n1"))


#' @rdname Design-class
#' @export
setGeneric("n2", function(d, x1, ...) standardGeneric("n2"))


#' @rdname Design-class
#' @export
setGeneric("n", function(d, x1, ...) standardGeneric("n"))


#' @rdname Design-class
#' @export
setGeneric("c2", function(d, x1, ...) standardGeneric("c2"))





#' Two-stage design
#'
#' \code{Design} is a abstract class for representing two-stage designs.
#'
#' @details Current main implementation: \code{\link{GQDesign-class}}
#' [TODO add some details on two-stage designs]
#'
#' @template DesignTemplate
#'
#' @exportClass Design
setClass("Design")


#' @describeIn Design must return numeric vector of length two giving early
#'     stopping for futility / efficacy boundaries; must be finite.
#'     Tries to access d@c1f and d@c1e otherwise throws an error.
#'     This means that the method need not be reimplemented as long as
#'     the subclass contains fields c1f and c2e.
#' @export
setMethod("early_stopping_bounds", signature("Design"),
    function(d, ...) {
        tryCatch(
            c(d@c1f, d@c1e),
            error = function(e) {
                stop("not implemented")
            }
        )
    })


#' @describeIn Design must return stage-one sample size of design d;
#'     default implementation tries to access d@n1 otherwise throws an error.
#' @export
setMethod("n1", signature("Design"),
          function(d, ...) {
              tryCatch(
                  d@n1,
                  error = function(e) {
                      stop("not implemented")
                  }
              )
          })


#' @describeIn Design stage-two sample size given stage one-outcome, must be implemented.
#' @export
setMethod("n2", signature("Design", "numeric"), function(d, x1, ...) stop("not implemented") )


#' @describeIn Design stage-two critical value given stage-one outcome, must be implemented.
#' @export
setMethod("c2", signature("Design", "numeric"), function(d, x1, ...) stop("not implemented") )


#' @describeIn Design overall sample size given stage-one outcome
#' @export
setMethod("n", signature("Design", "numeric"), function(d, x1, ...) n2(d, x1, ...) + n1(d, ...) )


#' @describeIn Design convert design object to numeric vector of design parameters, must be implemented.
#' @export
setMethod("as.numeric", signature("Design"), function(x) stop("not implemented"))


#' @describeIn Design convert vector of design parameters to design object of the
#'     same class as \code{object}; inverse to \code{as.numeric}, must be implemented.
#' @export
setMethod("update", signature("Design"), function(object) stop("not implemented"))
