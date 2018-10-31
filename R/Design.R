setClass("Design")

setGeneric("n1",
           function(d, ...) standardGeneric("n1")
)

setGeneric("n2",
        function(d, z1, ...) standardGeneric("n2")
)

setGeneric("c2",
        function(d, z1, ...) standardGeneric("c2")
)

setGeneric("conditional_power",
           function(d, z1, delta, ...) standardGeneric("conditional_power")
)

setGeneric("get_tunable_parameters", function(design, simplify = FALSE, ...) standardGeneric("get_tunable_parameters"))
