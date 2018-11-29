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
