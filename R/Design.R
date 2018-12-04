setClass("Design")


setGeneric("early_stopping_bounds", function(d, ...) standardGeneric("early_stopping_bounds"))
setMethod("early_stopping_bounds", signature("Design"),
    function(d, ...) {
        tryCatch(
            c(d@c1f, d@c1e),
            error = function(e) {
                stop("not implemented")
            }
        )
    })


setGeneric("n1", function(d, ...) standardGeneric("n1"))



setGeneric("n2", function(d, z1, ...) standardGeneric("n2"))



setGeneric("n", function(d, z1, ...) standardGeneric("n"))

setMethod("n", signature("Design", "numeric"),
    function(d, z1, ...) n2(d, z1, ...) + n1(d, ...) )



setGeneric("c2", function(d, z1, ...) standardGeneric("c2"))



setMethod("as.numeric", signature("Design"), function(x) stop("not implemented"))

setMethod("update", signature("Design"), function(object) stop("not implemented"))
