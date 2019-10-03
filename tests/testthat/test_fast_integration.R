context("fast integration for cdf/pdf")

test_that("normal priors and normal data build a correct data model", {

    library(dplyr)
    # TODO: Delete this when we do no longer depend on dplyr
    '
    # define scenarios
    mu  <- c(-.3, 0, .3, 1)
    sd  <- c(.1, .3, .5)
    n_1 <- c(10, 250, 900)
    # the entire grid takes way too much time!
    '
    mu  <- .3
    sd  <- .2
    n_1 <- 250

    for(i in 1:length(mu)) {
        for(j in 1:length(sd)) {
            for(k in 1:length(n_1)) {
                mu_prior <- mu[i]
                sd_prior <- sd[j]
                n        <- n_1[k]

                datadist <- Normal(two_armed = FALSE)
                prior    <- ContinuousPrior(
                    function(theta) dnorm(theta, mean = mu_prior, sd = sd_prior),
                    support = c(-5, 5),
                    tighten_support = T
                )


                m <- StageWiseDataModel(datadist, prior)


                # test marginal pdfs
                marginal_pdf_new  <- marginal_pdf(m, n = n, x = m@x)
                marginal_pdf_true <- dnorm(m@x, mean = mu_prior * sqrt(n), sd = sqrt(1 + n * sd_prior^2))
                # xbar | theta ~ N(theta, 1/n)
                # xbar         ~ N(mu_theta, (1 / n + sigma_theta^2))
                # z            ~ N(sqrt(n)* mu_theta, (1 + n * sigma_theta^2))

                testthat::expect_true(
                    all(abs(marginal_pdf_new - marginal_pdf_true) <= .025)
                )

                # test posterior pdfs
                xobs              <- sqrt(n) * 3 * mu_prior
                theta             <- seq(prior@support[1], prior@support[2], length.out = 1000)
                posterior_pdf_new <- posterior_pdf(m, theta = theta, n, x = xobs)

                mean_posterior     <- 1 / (1 / sd_prior^2 + n) * (mu_prior / sd_prior^2 + sqrt(n) * xobs)
                sd_posterior       <- sqrt(1 / (1 / sd_prior^2 + n))
                posterior_pdf_true <- dnorm(theta, mean = mean_posterior, sd = sd_posterior)

                testthat::expect_true(
                    all(abs(posterior_pdf_new - posterior_pdf_true) <= (.1 * max(posterior_pdf_true)))
                )

            }
        }
    }

})
