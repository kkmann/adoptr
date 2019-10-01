context("fast integration for cdf/pdf")

# p(x | n, theta)
Normal2 <- function(support_n, rel_support_x = c(-5, 5)) {

    res <- function(n, x, theta) {
        z <- pnorm(theta + rel_support_x[2], theta, 1/sqrt(n)) - pnorm(theta + rel_support_x[1], theta, 1/sqrt(n))
        ifelse(
            x <= theta + rel_support_x[2] & x >= theta + rel_support_x[1],
            dnorm(x, theta, 1/sqrt(n)) / z,
            0
        )
    }
    attr(res, "support_n") <- support_n
    attr(res, "rel_support_x") <- rel_support_x
    # attr(res, "class") <- c("pdf_stage", "function")

    return(res)
}

# p(theta)
NormalPrior2 <- function(mu, sd, support) {
    z   <- pnorm(support[2], mu, sd) - pnorm(support[1], mu, sd)
    res <- function(theta) dnorm(theta, mu, sd) / z
    attr(res, "support") <- support
    # attr(res, "class") <- c("pdf_theta", "function")

    return(res)
}

# new
pdf_x     <- Normal2(c(1, 100))
pdf_theta <- NormalPrior2(.3, .1, c(0, 1))
m         <- StageWiseDataModel(pdf_x, pdf_theta, dims = c(100, 500, 100))

new_marginal_pdf <- marginal_pdf(m, n = 10, x = m@x)

# old
datadist <- Normal(two_armed = FALSE)
prior    <- ContinuousPrior(pdf_theta, support = c(0, 1))

old_marginal_pdf <- predictive_pdf(datadist, prior, x1 = m@x, n1 = 10)


testthat::expect_true(
    all(abs(old_marginal_pdf - new_marginal_pdf) <= 10^(-4))
)


# SAVE FOR LATER
#
# posterior_pdf(m, theta = m@theta, 10, x = -1) %>%
#     {plot(m@theta, ., 'l')}
#
# x <- runif(100, -5, 10) %>% sort
# marginal_pdf(m, n = 10, x = x) %>%
#     {plot(x, ., 'l', xlim = c(-5, 10))}
#
# posterior_pdf(m, 0, 25, x = x) %>%
#     {plot(x, ., 'l', xlim = c(-2, 2))}
