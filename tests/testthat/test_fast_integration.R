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


sd_prior          <- .1
mu_prior          <- .3
support_prior     <- c(-3, 3)

n                 <- 25

# new
pdf_x             <- Normal2(c(1, 100))
pdf_theta         <- NormalPrior2(mu_prior, sd_prior, support_prior)
m                 <- StageWiseDataModel(pdf_x, pdf_theta, dims = c(100, 500, 100))

marginal_pdf_new  <- marginal_pdf(m, n = n, x = m@x)
marginal_pdf_true <- dnorm(m@x, mean = mu_prior * sqrt(n), sd = sqrt(1 + sd_prior^2))
# xbar | theta ~ N(theta, 1/n)
# xbar         ~ N(mu_theta, (1 + sigma_theta^2)/n)
# z            ~ N(sqrt(n)* mu_theta, (1 + sigma_theta^2))

testthat::expect_true(
    all(abs(marginal_pdf_new - marginal_pdf_true) <= .025)
)

# old
datadist <- Normal(two_armed = FALSE)
prior    <- ContinuousPrior(pdf_theta, support = support_prior)

marginal_pdf_old <- predictive_pdf(datadist, prior, x1 = m@x, n1 = n)

testthat::expect_true(
    all(abs(marginal_pdf_old - marginal_pdf_true) <= .025)
)

plot(m@x, marginal_pdf_true, 'l')
lines(m@x, marginal_pdf_new, col = 'red')
lines(m@x, marginal_pdf_old, col = 'blue')


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
