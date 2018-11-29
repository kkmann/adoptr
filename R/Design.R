setClass("Design")

setGeneric("n1",
           function(d, ...) standardGeneric("n1")
)

setGeneric("n2",
        function(d, z1, ...) standardGeneric("n2")
)

setGeneric("n",
           function(d, z1, ...) standardGeneric("n")
)

setMethod("n", signature("Design", "numeric"),
          function(d, z1, ...) n2(d, z1, ...) + n1(d, ...)
)

setGeneric("c2",
        function(d, z1, ...) standardGeneric("c2")
)

setGeneric("conditional_power",
           function(d, z1, delta, ...) standardGeneric("conditional_power")
)
