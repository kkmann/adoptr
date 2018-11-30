setClass("ConditionalScore", representation(prior = "Distribution"))

# evaluate a score
setGeneric("eval", function(s, design, ...) standardGeneric("eval"))
# allow cusotm implementation depending on design
setGeneric(".eval_specific", function(s, design, ...) standardGeneric(".eval_specific"))

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
setMethod("eval", signature("IntegralScore", "Design"),
    function(s, design, specific = TRUE, ...) {
        # TODO: currently ignores the possibility of early stopping/uncontinuus
        # conditional scores - might get better when checking for early stopping
        # and integrating separately!
        if (specific) { # use design-specific implementation
            return(.eval_specific(s, design, ...))
        } else { # use generic approach
            # integrand is the conditional score as function of z1 times the
            # predictive pdf given the scores prior
            integrand <- function(z1) eval(s@conditional_score, design, z1, ...) *
                predictive_pdf(s@conditional_score@prior, z1, n1(design), ...)
            z1_bounds <- get_z1_bounds(s@conditional_score@prior, mass = .99998)
            # use adaptive quadrature to integrate - only relies on generic interface
            # provided by 'Design', no special optimization for particular
            # design implementation
            return(stats::integrate(integrand, z1_bounds[1], z1_bounds[2])$value)
        }
    })

setMethod(".eval_specific", signature("IntegralScore", "Design"),
    function(s, design, ...) stop("not implemented") )
