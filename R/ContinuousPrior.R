# setClass("Distribution")
#
# setGeneric("condition", function(dist, lower, upper, ...) standardGeneric("condition"))
#
# setGeneric("posterior", function(dist, z1, n1, ...) standardGeneric("posterior"))
#
# setGeneric("predictive_pdf", function(dist, z1, n1, ...) standardGeneric("predictive_pdf"))
#
# setGeneric("predictive_cdf", function(dist, z1, n1, ...) standardGeneric("predictive_cdf"))
#
# setGeneric("expectation", function(dist, f, ...) standardGeneric("expectation"))
#
# setGeneric("get_z1_bounds", function(dist, ...) standardGeneric("get_z1_bounds"))
#
#
#
# setClass("PointMassDistribution",
#          representation(
#                  theta = "numeric",
#                  mass  = "numeric"
#          ),
#          contains = "Distribution"
# )
#
#
#
# PointMassDistribution <- function(theta, mass) {
#         if (sum(mass) != 1)
#                 stop("mass must sum to one")
#         new("PointMassDistribution", theta = theta, mass = mass)
# }
#
# setMethod("get_z1_bounds",
#           signature("PointMassDistribution"),
#           function(dist, lower, mass = .99998, ...) {
#                   return(c(
#                         min(dist@theta) + qnorm((1 - mass) / 2),
#                         max(dist@theta) - qnorm((1 - mass) / 2)
#                   ))
#           }
# ) # condition
#
# setMethod("condition",
#           signature("PointMassDistribution", "numeric", "numeric"),
#           function(dist, lower, upper, ...) {
#                   idx  <- (dist@theta >= lower) & (dist@theta <= upper)
#                   return(PointMassDistribution(
#                           dist@theta[idx], dist@mass[idx]/sum(dist@mass[idx])
#                   ))
#           }
# ) # condition
#
# setMethod("predictive_pdf",
#           signature("PointMassDistribution", "numeric"),
#           function(dist, z1, n1, ...) {
#                   k   <- length(dist@theta)
#                   res <- numeric(length(z1))
#                   for (i in 1:k) {
#                           res <- res + dist@mass[i] * dnorm(z1, mean = sqrt(n1) * dist@theta, sd = 1)
#                   }
#                   return(res)
#           }
# ) # predictive_pdf
#
#
# setMethod("predictive_cdf",
#           signature("PointMassDistribution", "numeric"),
#           function(dist, z1, n1, ...) {
#                   k   <- length(dist@theta)
#                   res <- numeric(length(z1))
#                   for (i in 1:k) {
#                           res <- res + dist@mass[i] * pnorm(z1, mean = sqrt(n1) * dist@theta, sd = 1)
#                   }
#                   return(res)
#           }
# ) # predictive_cdf
#
#
# setMethod("posterior",
#           signature("PointMassDistribution", "numeric"),
#           function(dist, z1, n1, ...) {
#                   mass <- dist@mass * dnorm(z1, mean = sqrt(n1) * dist@theta, sd = 1)
#                   mass <- mass / sum(mass) # normalize
#                   return(PointMassDistribution(dist@theta, mass))
#           }
# ) # posterior
#
# setMethod("expectation",
#           signature("PointMassDistribution", "function"),
#           function(dist, f, ...) {
#                   ff <- sapply(dist@theta, f, ...)
#                   return(sum(dist@mass * ff))
#           }
# ) # expectation
