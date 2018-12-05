setClass("ConditionalPower", contains = "ConditionalScore")

ConditionalPower <- function(dist, prior) new("ConditionalPower", distribution = dist, prior = prior)

setMethod("evaluate", signature("ConditionalPower", "Design"),
          function(s, design, x1, ...) {
              sapply(x1,
                  function(x1) expectation(
                      posterior(s@distribution, s@prior, x1, n1(design), ...),
                      function(theta)
                          1 - cumulative_distribution_function(s@distribution, c2(design, x1), n2(design, x1), theta)
                  )
              )
          })
