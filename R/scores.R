setClass("ESS", representation(prior = "Prior"), contains = "Score")

ESS <- function(prior) {
        new("ESS", prior = prior)
}

setMethod("eval", signature("ESS", "Design"), function(s, design, ...) {
        mass <- cubature::hcubature(
                f = function(x) dnorm(x[1], mean = x[2], sd = 1) * s@prior@p(x[2]),
                lowerLimit = c(-1, -1),
                upperLimit = c(4, 1),
                absError   = 1e-3
        )
        res <- cubature::hcubature(
                f = function(x) n2(design, x[1]) * dnorm(x[1], mean = x[2], sd = 1) * s@prior@p(x[2]),
                lowerLimit = c(-1, -1),
                upperLimit = c(4, 1),
                absError   = .1
        )$integral / mass$integral
        return(res + n1(design))
})

setClass("ESS2", representation(delta = "numeric"), contains = "Score")

ESS2 <- function(delta) {
        new("ESS2", delta = delta)
}

setMethod("eval", signature("ESS2", "Design"), function(s, design, ...) {
        xx   <- seq(-1, 4, by = .001)
        mass <- .001 * sum(dnorm(xx, mean = s@delta*sqrt(n1(design)), sd = 1))
        res  <- .001 * sum(n2(design, xx) * dnorm(xx, mean = s@delta*sqrt(n1(design)), sd = 1)) / mass
        return(res + n1(design))
})

# setClass("ExpectedPower", representation(prior = "Prior"), contains = "Score")
#
# ExpectedPower <- function(prior) {
#         new("ExpectedPower", prior = prior)
# }
#
# setMethod("eval", signature("ExpectedPower", "Design"), function(s, design, ...) {
#         mass <- cubature::hcubature(
#                 f = function(x) dnorm(x[1], mean = x[2], sd = 1) * s@prior@p(x[2]),
#                 lowerLimit = c(-1, -1),
#                 upperLimit = c(4, 1),
#                 absError   = 1e-3
#         )
#         res <- cubature::hcubature(
#                 f = function(x) (1 - pnorm(c2(design, x[1]) - sqrt(n2(design, x[1])))) * dnorm(x[1], mean = x[2], sd = 1) * s@prior@p(x[2]),
#                 lowerLimit = c(-1, -1),
#                 upperLimit = c(4, 1),
#                 absError   = .1
#         )$integral / mass$integral
#         return(res)
# })


setClass("Power", representation(delta = "numeric"), contains = "Score")

Power <- function(delta) {
        new("Power", delta = delta)
}

setMethod("eval", signature("Power", "Design"), function(s, design, ...) {
        xx   <- seq(-1, 4, by = .001)
        mass <- .001 * sum(dnorm(xx, mean = s@delta*sqrt(n1(design)), sd = 1))
        # res <- integrate(
        #         function(z1) {
        #                 (1 - pnorm(c2(design, z1) - sqrt(n2(design, z1))*s@delta)) * dnorm(z1, mean = s@delta*sqrt(n1(design)), sd = 1)
        #         },
        #         -1,
        #         4
        # )
        res  <- .001 * sum(
                (1 - pnorm(c2(design, xx) - sqrt(n2(design, xx))*s@delta)) * dnorm(xx, mean = s@delta*sqrt(n1(design)), sd = 1)
        ) / mass
        return(res)
})
