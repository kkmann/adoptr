setClass("ConditionalScore", representation(prior = "Distribution"))

# to be defined by each score, eval(s, design, z1, ...) for conditional scores
setGeneric("eval",
           function(s, design, ...) standardGeneric("eval")
)

# better name? integrate a conditional score with respect to Z1
setGeneric("integrate",
           function(s, ...) standardGeneric("integrate")
)

# essentialy just changes class and stores conditional score to allow
# different implementation of 'eval'
setMethod("integrate",
          signature("ConditionalScore"),
          function(s, ...) {
                  new("IntegralScore", conditional_score = s)
          }
)

# well, result of 'integrate'
setClass("IntegralScore", representation(conditional_score = "ConditionalScore"))

# generic evaluation of integral score, this is where the problems lie ;)
# uses stats::integrate to integrate conditional score over Z1, not working for
# optimization - need custom implementation for signature("IntegralScore", "BSDesign")
setMethod("eval",
          signature("IntegralScore", "Design"),
          function(s, design, ...) {
                  integrand <- function(z1) eval(s@conditional_score, design, z1, ...) *
                          predictive_pdf(s@conditional_score@prior, z1, n1(design), ...)
                  z1_bounds <- get_z1_bounds(s@conditional_score@prior, mass = .99998)
                  return(stats::integrate(integrand, z1_bounds[1], z1_bounds[2])$value)
          }
)
