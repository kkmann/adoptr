setClass("Problem", representation(minimize = "Score", subject_to = "list"))

Problem <- function(minimize, subject_to) {
        new("Problem", minimize = minimize, subject_to = subject_to)
}

setGeneric("get_Lagrangian", function(problem) standardGeneric("get_Lagrangian"))

setMethod("get_Lagrangian", signature("Problem"), function(problem) {
        function(design, multipliers) {
                res <- eval(problem@minimize, design)
                for (i in 1:length(problem@subject_to)) {
                        cnsrt   <- problem@subject_to[[i]]
                        penalty <- multipliers[i] * (eval(cnsrt@score, design) - cnsrt@RHS) # score <= rhs
                        if (cnsrt@orientation == ">=")
                                penalty <- penalty * (-1)
                        res   <- res + penalty
                }
                return(res)
        }
})
