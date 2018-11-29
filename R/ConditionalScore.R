setClass("ConditionalScore", representation(prior = "Distribution"))

setGeneric("integrate",
           function(s, ...) standardGeneric("integrate")
)

setClass("IntegralScore", representation(conditional_score = "ConditionalScore"), contains = "Score")

setMethod("integrate",
          signature("ConditionalScore"),
          function(s, ...) {
                new("IntegralScore", conditional_score = s)
          }
)

setMethod("eval",
          signature("IntegralScore", "Design"),
          function(s, design, ...) {
                integrand <- function(z1) eval(s@conditional_score, design, z1, ...) *
                        predictive_pdf(s@conditional_score@prior, z1, n1(design), ...)
                z1_bounds <- get_z1_bounds(s@conditional_score@prior, mass = .99998)
                return(stats::integrate(integrand, z1_bounds[1], z1_bounds[2])$value)
          }
)
