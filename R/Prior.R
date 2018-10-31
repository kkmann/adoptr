setClass("Prior", representation(p = "function"))

Prior <- function(p) {
        mass <- integrate(p, -3, 3)$value
        if (mass < .99) warning(sprintf("prior mass in (-3, 3) is low (%.2f)", mass))
        new("Prior", p = p)
}

setGeneric("condition", function(prior, lower, upper, ...) standardGeneric("condition"))

setMethod("condition",
          signature("Prior", "numeric", "numeric"),
          function(prior, lower, upper, abs.error.tol = 1e-4, ...) {
                  mass <- integrate(prior@p, lower, upper)
                  if (mass$abs.error > abs.error.tol)
                          warning(sprintf("absolute error of %.2e exceeds tolerance %.2e", mass$abs.error, abs.error.tol))
                  p <- function(delta) ifelse(delta >= lower & delta < upper, 1, 0) * prior@p(delta) / mass$value
                  new("Prior", p = p)
          }
) # condition

setGeneric("prior_predictive", function(prior, z1, ...) standardGeneric("prior_predictive"))

setMethod("prior_predictive",
          signature("Prior", "numeric"),
          function(prior, z1, rel.error.tol = 1e-2, lower = -Inf, upper = Inf, ...) {
                grid <- expand.grid(z1 = z1, delta = seq(-3, 3, .01))
                unnormalized <- rowSums(matrix(
                        apply(grid, 1, function(x) dnorm(x[1], mean = x[2]) * prior@p(x[2])),
                        nrow = length(z1)
                ))
                return(unnormalized / sum(unnormalized))
          }
)
