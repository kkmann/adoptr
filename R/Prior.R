setClass("Prior")

# all priors should be limited, return bounds
setGeneric("bounds", function(prior, ...) standardGeneric("bounds"))

# expected value of function w.r.t. prior
setGeneric("expectation", function(prior, f, ...) standardGeneric("expectation"))

# conditionl prior on interval c(lower, upper)
setGeneric("condition", function(prior, interval, ...) standardGeneric("condition"))

# predictive distribution (assuming normal model)
setGeneric("predictive_pdf", function(prior, z1, n1, ...) standardGeneric("predictive_pdf"))

setGeneric("predictive_cdf", function(prior, z1, n1, ...) standardGeneric("predictive_cdf"))

# posterior after observing sufficient statistic (z, n)
setGeneric("posterior", function(prior, z1, n1, ...) standardGeneric("posterior"))


