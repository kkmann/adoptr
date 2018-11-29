setClass("BSDesign",
        representation(
                n1         = "numeric",
                early_stopping_for_futility = "numeric",
                early_stopping_for_efficacy = "numeric",
                weights_n2 = "numeric",
                weights_c2 = "numeric",
                n_knots = "integer"
        ),
        contains = "Design"
)

BSDesign <- function(
        n1,
        early_stopping_for_futility,
        early_stopping_for_efficacy,
        weights_n2,
        weights_c2,
        n_knots
) {

        new("BSDesign",
                n1 = n1,
                early_stopping_for_futility = early_stopping_for_futility,
                early_stopping_for_efficacy = early_stopping_for_efficacy,
                weights_n2 = weights_n2,
                weights_c2 = weights_c2,
                n_knots = n_knots
        )

}


setMethod("update",
          signature("BSDesign"),
          function(object, tunable_parameters, ...) {
                  k <- length(tunable_parameters)
                  if ((k - 3) %% 2 != 0)
                          stop("parameter length does not fit")
                  n_weights <- (k - 3)/2
                  new("BSDesign",
                      n1 = tunable_parameters[1],
                      early_stopping_for_futility = tunable_parameters[2],
                      early_stopping_for_efficacy = tunable_parameters[3],
                      weights_n2 = tunable_parameters[4:(3 + n_weights)],
                      weights_c2 = tunable_parameters[(4 + n_weights):length(tunable_parameters)],
                      n_knots = object@n_knots
                  )
          }
)

setMethod("n1",
          signature("BSDesign"),
          function(d, ...) d@n1
)


setGeneric("get_knots",
           function(d, ...) standardGeneric("get_knots")
)

setMethod("get_knots",
          signature("BSDesign"),
          function(d, ...) seq(d@early_stopping_for_futility, d@early_stopping_for_efficacy, length.out = d@n_knots) # TODO: make sure > 0
)

setMethod("n2",
        signature("BSDesign", "numeric"),
        function(d, z1, ...) {
                return(
                        ifelse(z1 < d@early_stopping_for_futility | z1 > d@early_stopping_for_efficacy, 0, 1) *
                                pmax(0, spline(get_knots(d), d@weights_n2, xout = z1, method = "natural")$y)
                )
        }
)

setMethod("c2",
        signature("BSDesign", "numeric"),
        function(d, z1, ...) {
                return(
                        spline(get_knots(d), d@weights_c2, xout = z1, method = "natural")$y *
                        ifelse(z1 < d@early_stopping_for_futility, Inf, 1) *
                        ifelse(z1 > d@early_stopping_for_efficacy, -Inf, 1)
                )
        }
)

setMethod("conditional_power",
          signature("BSDesign", "numeric", "numeric"),
          function(d, z1, delta, ...) {
                  n2 <- n2(d, z1)
                  c2 <- c2(d, z1)
                  return(1 - pnorm(c2 - sqrt(n2)*delta))
          }
)

setMethod("as.numeric",
          signature("BSDesign"),
          function(x, ...) {
                  return(c(
                          x@n1,
                          x@early_stopping_for_futility,
                          x@early_stopping_for_efficacy,
                          x@weights_n2,
                          x@weights_c2
                        )
                  )
          }
)

# need to do this differently: piecewiese Gauss-Quadrature with 3 pivots ->
# # exact for piecewise polynomials of degree 5, 3 sections -> 9 pivots
# setMethod("eval",
#           signature("IntegralScore", "BSDesign"),
#           function(s, design, ...) {
#                   integrand <- function(z1) eval(s@conditional_score, design, z1, ...) *
#                           predictive_pdf(s@conditional_score@prior, z1, n1(design), ...)
#                   get_knots(d)
#                   z1_bounds <- get_z1_bounds(s@conditional_score@prior, mass = .99998)
#                   return(stats::integrate(integrand, z1_bounds[1], z1_bounds[2])$value)
#           }
# )

