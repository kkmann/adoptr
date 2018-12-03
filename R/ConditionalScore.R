setClass("ConditionalScore", representation(prior = "Prior"))

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
