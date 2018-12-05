setClass("Design")

#' @export
setGeneric("early_stopping_bounds", function(d, ...) standardGeneric("early_stopping_bounds"))
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

#' @export
setGeneric("n1", function(d, ...) standardGeneric("n1"))


#' @export
setGeneric("n2", function(d, x1, ...) standardGeneric("n2"))


#' @export
setGeneric("n", function(d, x1, ...) standardGeneric("n"))

setMethod("n", signature("Design", "numeric"),
    function(d, x1, ...) n2(d, x1, ...) + n1(d, ...) )

#' @export
setGeneric("c2", function(d, x1, ...) standardGeneric("c2"))

#' @export
setGeneric("get_knots", function(d, ...) standardGeneric("get_knots"))

#' @export
setMethod("as.numeric", signature("Design"), function(x) stop("not implemented"))

#' @export
setMethod("update", signature("Design"), function(object) stop("not implemented"))
