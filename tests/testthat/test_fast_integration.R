context("fast integration for cdf/pdf")

# p(x | n, theta)
Normal2 <- function(support_n, rel_support_x = c(-5, 5)) {

    res <- function(n, x, theta) {
        z <- pnorm(sqrt(n) * theta + rel_support_x[2], sqrt(n) * theta, 1) - pnorm(sqrt(n) * theta + rel_support_x[1], sqrt(n) * theta, 1)

        ifelse(
            x <= sqrt(n) * theta + rel_support_x[2] & x >= sqrt(n) * theta + rel_support_x[1],
            dnorm(x, sqrt(n) * theta, 1) / z,
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

n                 <- 200 # play around with that here

# new
data_pdf          <- Normal2(c(1, 250))
prior_pdf         <- NormalPrior2(mu_prior, sd_prior, support_prior)
m                 <- StageWiseDataModel(data_pdf, prior_pdf)

marginal_pdf_new  <- marginal_pdf(m, n = n, x = m@x)
marginal_pdf_true <- dnorm(m@x, mean = mu_prior * sqrt(n), sd = sqrt(1 + n * sd_prior^2))
# xbar | theta ~ N(theta, 1/n)
# xbar         ~ N(mu_theta, (1 + sigma_theta^2)/n)
# z            ~ N(sqrt(n)* mu_theta, (1 + sigma_theta^2))

testthat::expect_true(
    all(abs(marginal_pdf_new - marginal_pdf_true) <= .025)
)

# old
datadist <- Normal(two_armed = FALSE)
prior    <- ContinuousPrior(prior_pdf, support = support_prior, tighten_support = T)

marginal_pdf_old <- predictive_pdf(datadist, prior, x1 = m@x, n1 = n)

testthat::expect_true(
    all(abs(marginal_pdf_old - marginal_pdf_true) <= .025)
)

plot(m@x, marginal_pdf_true, 'l')
lines(m@x, marginal_pdf_new, col = 'red')
lines(m@x, marginal_pdf_old, col = 'blue')


# posteriors
xobs              <- sqrt(n) * 3 * mu_prior
theta             <- seq(-3, 3, length.out = 1000)
posterior_pdf_new <- posterior_pdf(m, theta = theta, n, x = xobs)
posterior_pdf_old <- posterior(datadist, prior, x1 = xobs, n1 = n)@pdf(theta)

mean_posterior     <- 1 / (1 / sd_prior^2 + n) * (mu_prior / sd_prior^2 + sqrt(n) * xobs)
sd_posterior       <- sqrt(1 / (1 / sd_prior^2 + n))
posterior_pdf_true <- dnorm(theta, mean = mean_posterior, sd = sd_posterior)


plot(theta, posterior_pdf_true, 'l', ylim = c(-1, 10))
lines(theta, posterior_pdf_new, col = 'red')
lines(theta, posterior_pdf_old, col = 'blue')
