setClass("FivePointDesign",
        representation(
                n1         = "numeric",
                early_stopping_for_futility = "numeric",
                early_stopping_for_efficacy = "numeric",
                n2_pivots = "numeric",
                c2_pivots = "numeric"
        ),
        contains = "Design"
)

FivePointDesign <- function(
        n1,
        early_stopping_for_futility,
        early_stopping_for_efficacy,
        n2_pivots,
        c2_pivots
) {
        if (length(n2_pivots) != 5 | length(c2_pivots) != 5)
                stop("both pivot vectors must have length 5")
        new("FivePointDesign",
                n1 = n1,
                early_stopping_for_futility = early_stopping_for_futility,
                early_stopping_for_efficacy = early_stopping_for_efficacy,
                n2_pivots = n2_pivots,
                c2_pivots = c2_pivots
        )

}


setMethod("update",
          signature("FivePointDesign"),
          function(object, tunable_parameters, ...) {
                  k <- length(tunable_parameters)
                  if (k != 13)
                          stop("parameter length does not fit")
                  n_pivots <- 5
                  new("FivePointDesign",
                      n1 = tunable_parameters[1],
                      early_stopping_for_futility = tunable_parameters[2],
                      early_stopping_for_efficacy = tunable_parameters[3],
                      n2_pivots = tunable_parameters[4:(3 + n_pivots)],
                      c2_pivots = tunable_parameters[(4 + n_pivots):length(tunable_parameters)]
                  )
          }
)

setMethod("n1",
          signature("FivePointDesign"),
          function(d, ...) d@n1
)


setMethod("get_knots",
          signature("FivePointDesign"),
          function(d, ...) seq(d@early_stopping_for_futility, d@early_stopping_for_efficacy, length.out = 5)
)

setMethod("n2",
        signature("FivePointDesign", "numeric"),
        function(d, z1, ...) {
                return(
                        ifelse(z1 < d@early_stopping_for_futility | z1 > d@early_stopping_for_efficacy, 0, 1) *
                                pmax(0, approx(get_knots(d), d@n2_pivots, xout = z1, method = "linear", rule = 2)$y)
                )
        }
)

setMethod("c2",
        signature("FivePointDesign", "numeric"),
        function(d, z1, ...) {
                return(
                        approx(get_knots(d), d@c2_pivots, xout = z1, method = "linear", rule = 2)$y *
                        ifelse(z1 < d@early_stopping_for_futility, Inf, 1) *
                        ifelse(z1 > d@early_stopping_for_efficacy, -Inf, 1)
                )
        }
)

setMethod("conditional_power",
          signature("FivePointDesign", "numeric", "numeric"),
          function(d, z1, delta, ...) {
                  n2 <- n2(d, z1)
                  c2 <- c2(d, z1)
                  return(1 - pnorm(c2 - sqrt(n2)*delta))
          }
)

setMethod("as.numeric",
          signature("FivePointDesign"),
          function(x, ...) {
                  return(c(
                          x@n1,
                          x@early_stopping_for_futility,
                          x@early_stopping_for_efficacy,
                          x@n2_pivots,
                          x@c2_pivots
                        )
                  )
          }
)


# need to do this differently: piecewiese Gauss-Quadrature with 3 pivots ->
# # exact for piecewise polynomials of degree 5, 3 sections -> 9 pivots
setMethod("eval",
          signature("IntegralScore", "FivePointDesign"),
          function(s, design, ...) {
                  # early stopping probabilities (conditional score constant here!:
                  poef <- predictive_cdf(s@conditional_score@prior,
                                         design@early_stopping_for_futility, n1(design))
                  poee <- 1 - predictive_cdf(s@conditional_score@prior,
                                         design@early_stopping_for_efficacy, n1(design))
                  # continuation region
                  integrand   <- function(z1) eval(s@conditional_score, design, z1, ...) *
                          predictive_pdf(s@conditional_score@prior, z1, n1(design), ...)
                  weights     <- c(7, 32, 12, 32, 7)
                  h           <- (design@early_stopping_for_efficacy - design@early_stopping_for_futility)/4
                  mid_section <- 2/45 * h * sum(weights * integrand(get_knots(design)))
                  # compose
                  res <- poef * eval(
                                s@conditional_score, design,
                                design@early_stopping_for_futility - sqrt(.Machine$double.eps) # slightly smaller than stopping for futility
                          ) +
                          mid_section +
                          poee * eval(
                                  s@conditional_score, design,
                                  design@early_stopping_for_efficacy + sqrt(.Machine$double.eps)
                          )
                 return(res)
          }
)

