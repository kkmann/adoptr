setClass("ConditionalPower", contains = "ConditionalScore")

ConditionalPower <- function(prior) new("ConditionalPower", prior = prior)

setMethod("evaluate", signature("ConditionalPower", "Design"),
          function(s, design, x1, ...) {
                res <- sapply(x1, function(x1) expectation(
                                posterior(s@prior, x1, n1(design), ...),
                                function(theta) conditional_power(design, x1, theta)
                        )
                )
                return(res)
          })
