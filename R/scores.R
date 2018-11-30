setClass("SampleSize", contains = "ConditionalScore")

SampleSize <- function(prior) {
        new("SampleSize", prior = prior)
}

setMethod("eval",
          signature("SampleSize", "Design"),
          function(s, design, z1, ...) n2(design, z1) + n1(design)
)




setClass("ConditionalPower", contains = "ConditionalScore")

ConditionalPower <- function(prior) {
        new("ConditionalPower", prior = prior)
}

setMethod("eval",
          signature("ConditionalPower", "Design"),
          function(s, design, z1, ...) {
                res <- sapply(z1, function(z1) expectation(
                                posterior(s@prior, z1, n1(design), ...),
                                function(theta) conditional_power(design, z1, theta)
                        )
                )
                return(res)
          }
)


# smoothness terms must be implemented 'by-hand' since they are independent of any
# distributions!
setClass("Smoothness_n2", representation(
        h = "numeric"
    ),
    contains = "IntegralScore")

Smoothness_n2 <- function(prior, h = sqrt(.Machine$double.eps)) {
    new("Smoothness_n2", h = h)
}

setMethod("eval", signature("Smoothness_n2", "Design"),
          function(s, design, specific = TRUE, ...) {
              if (specific) { # use design-specific implementation
                  return(.eval_specific(s, design, ...))
              } else { # use generic approach
                  # integrand is the finite difference approximation of the
                  # squared derivative
                  integrand <- function(z1) (n2(design, z1 + s@h/2) - n2(design, z1 - s@h/2))^2
                  z1_bounds <- tryCatch(
                      early_stopping_bounds(design, ...),
                      error = function(e) { # assume no early stopping, integrate over the entire range!
                          qnorm(.0005, .9995, mean = bounds(s@conditional_score@prior) * n1(design), sd = 1)
                      }
                  )
                  # use adaptive quadrature to integrate - only relies on generic interface
                  # provided by 'Design', no special optimization for particular
                  # design implementation
                  return(1 / diff(z1_bounds) * stats::integrate(integrand, z1_bounds[1], z1_bounds[2])$value)
              }
          })
