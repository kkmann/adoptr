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



setGeneric("conditional_power", function(d, z1, delta, ...) standardGeneric("conditional_power"))

# TODO: allow null other than delta == 0!
setMethod("conditional_power", signature("Design", "numeric", "numeric"),
    function(d, z1, delta, ...) 1 - pnorm(c2(d, z1) - sqrt(n2(d, z1))*delta) )



setMethod("as.numeric", signature("Design"), function(x) stop("not implemented"))

setMethod("update", signature("Design"), function(object) stop("not implemented"))
