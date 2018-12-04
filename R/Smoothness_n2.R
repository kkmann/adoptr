# smoothness terms must be implemented 'by-hand' since they are independent of any
# distributions! (cannot be handeld by 'integrate')
# TODO, this can probably be solved in a nicer way ;)
setClass("Smoothness_n2", representation(
    h = "numeric"
    ),
    contains = "UnconditionalScore")

Smoothness_n2 <- function(h = sqrt(.Machine$double.eps)) new("Smoothness_n2", h = h)

setMethod("evaluate", signature("Smoothness_n2", "Design"),
          function(s, design, specific = TRUE, ...) {
              if (specific) { # use design-specific implementation
                  return(.evaluate(s, design, ...))
              } else {
                  # use generic approach
                  # integrand is the finite difference approximation of the
                  # squared derivative
                  integrand <- function(x1) ((n2(design, x1 + s@h/2) - n2(design, x1 - s@h/2))/s@h)^2
                  x1_bounds <- tryCatch( # only evaluate n2 within continuation region!
                      early_stopping_bounds(design, ...) + c(s@h, -s@h),
                      error = function(e) { # assume no early stopping, integrate over the entire range!
                          qnorm(c(.0005, .9995), mean = bounds(s@cs@prior) * n1(design), sd = 1)
                      }
                  )
                  # use adaptive quadrature to integrate - only relies on generic interface
                  # provided by 'Design', no special optimization for particular
                  # design implementation
                  return(1 / diff(x1_bounds) * stats::integrate(integrand, x1_bounds[1], x1_bounds[2])$value)
              }
          })
