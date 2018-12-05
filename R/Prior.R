setClass("Prior")

# all priors should be limited, return bounds
setGeneric("bounds", function(dist, ...) standardGeneric("bounds"))

# expected value of function w.r.t. prior
setGeneric("expectation", function(dist, f, ...) standardGeneric("expectation"))

# conditionl prior on interval c(lower, upper)
setGeneric("condition", function(dist, interval, ...) standardGeneric("condition"))

# predictive distribution (assuming normal model)
setGeneric("predictive_pdf", function(dist, prior, x1, n1, ...) standardGeneric("predictive_pdf"))

setGeneric("predictive_cdf", function(dist, prior, x1, n1, ...) standardGeneric("predictive_cdf"))

# posterior after observing sufficient statistic (z, n)
setGeneric("posterior", function(dist, prior, x1, n1, ...) standardGeneric("posterior"))
