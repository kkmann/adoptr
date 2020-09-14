context("Binomial data distribution")


test_that("Constructor works", {

    expect_true(
        Binomial(.2)@two_armed)

    expect_true(
        !Binomial(.4, two_armed = FALSE)@two_armed)

    expect_error(
        Binomial(-.1))

    expect_error(
        Binomial(1.1))

}) # end 'constructor works'


test_that("pdf is defined correctly", {

    dist    <- Binomial(.1, TRUE)
    x       <- seq(-3, 3, by = .1)
    n       <- 22
    theta   <- .35
    r_A     <- theta + dist@rate_control
    sigma_A <- sqrt(sum(c(r_A * (1 - r_A), dist@rate_control * (1 - dist@rate_control))))
    sigma_0 <- sqrt(2 * mean(c(r_A, dist@rate_control)) * (1 - mean(c(r_A, dist@rate_control))))
    expect_equal(
        probability_density_function(dist, x, n, theta),
        stats::dnorm(x, mean = sqrt(n) * theta / sigma_0, sd = sigma_A / sigma_0),
        tolerance = 1e-6, scale = 1)
})



test_that("cdf is defined correctly", {

    dist  <- Binomial(.3, TRUE)
    x     <- seq(-3, 3, by = .1)
    n     <- 120
    theta <- .02
    r_A     <- theta + dist@rate_control
    sigma_A <- sqrt(sum(c(r_A * (1 - r_A), dist@rate_control * (1 - dist@rate_control))))
    sigma_0 <- sqrt(2 * mean(c(r_A, dist@rate_control)) * (1 - mean(c(r_A, dist@rate_control))))
    expect_equal(
        cumulative_distribution_function(dist, x, n, theta),
        stats::pnorm(x, mean = sqrt(n) * theta / sigma_0, sd = sigma_A / sigma_0),
        tolerance = 1e-6, scale = 1)
})



test_that("quantile is defined correctly", {

    dist  <- Binomial(.6, FALSE)
    probs <- seq(.01, 1, by = .01)
    n     <- 99
    theta <- .125
    r_A     <- theta + dist@rate_control
    sigma_A <- sqrt(r_A * (1 - r_A))
    sigma_0 <- sqrt(2 * r_A * (1 - r_A))
    expect_equal(
        quantile(dist, probs, n, theta),
        stats::qnorm(probs, mean = sqrt(n) * theta / sigma_0, sd = sigma_A / sigma_0),
        tolerance = 1e-6, scale = 1)

    expect_warning(
        quantile(dist, -1, n, theta)
    )

})


test_that("vectorization works", {
    dist  <- Binomial(.3)
    prior <- PointMassPrior(c(.3, .4), c(.4, .6))
    expect_equal(
        predictive_pdf(dist, prior, c(0, 1), 10),
        sapply(c(0, 1), function(x) predictive_pdf(dist, prior, x, 10))
    )
})


# Let us check if the statistical tests are working correctly
r_c    <- .4
r_e    <- .6
r_diff <- r_e - r_c
dist   <- Binomial(r_c)
alpha  <- 0.025
beta   <- 0.2
n      <- ceiling(rpact::getSampleSizeRates(alpha = alpha, beta = beta,
                                            groups = 2, normalApproximation = TRUE,
                                            pi1 = r_e, pi2 = r_c)$nFixed1)
c      <- qnorm(1 - alpha)



test_that("Approximation works in cdf", {
    iters   <- 1e6
    p_c     <- stats::rbinom(iters, n, r_c) / n
    p_e     <- stats::rbinom(iters, n, r_e) / n
    p_e_2   <- stats::rbinom(iters, n, r_c) / n
    sigma_A <- sqrt(p_e * (1 - p_e) + p_c * (1 - p_c))
    sigma_0 <- sqrt(p_e_2 * (1 - p_e_2) + p_c * (1 - p_c))

    expect_equal(
        mean(1 - cumulative_distribution_function(dist, sqrt(n) * (p_e - p_c) / sigma_A, n, 0) <= alpha),
        0.807, # reference: M. Kieser (2018), Fallzahlberechnung in der medizinischen Forschung, Tab. 5.2
        tolerance = 1e-2, scale = 1
    )

    expect_equal(
        mean(1 - cumulative_distribution_function(dist, sqrt(n) * (p_e_2 - p_c) / sigma_0, n, 0) <= alpha),
        alpha,
        tolerance = 1 / sqrt(iters), scale = 1
    )
})


test_that("Approximation works in evaluate", {
    design <- OneStageDesign(n, c)
    H_0    <- PointMassPrior(0, 1)
    H_1    <- PointMassPrior(r_diff, 1)

    expect_equal(
        evaluate(Power(dist, H_1), design),
        1 - beta,
        tolerance = 1e-3, scale = 1
    )

    expect_equal(
        evaluate(Power(dist, H_0), design),
        alpha,
        tolerance = 1e-3, scale = 1
    )

})




test_that("Approximation works in simulate", {
    iters <- 1e6
    test_statistic_H0 <- simulate(object = dist, nsim = iters, n = n, theta = 0, seed = 42)
    test_statistic_H1 <- simulate(object = dist, nsim = iters, n = n, theta = r_diff, seed = 42)

    expect_equal(
        mean(test_statistic_H1 > c),
        1 - beta,
        tolerance = 1 / sqrt(iters), scale = 1
    )

    expect_equal(
        mean(test_statistic_H0 > c),
        alpha,
        tolerance = 1 / sqrt(iters), scale = 1
    )

})



test_that("show method", {

    expect_equal(
        capture.output(show(Binomial(.2))),
        "Binomial<two-armed>, response rate in control group: 0.2 "
    )

    expect_equal(
        capture.output(show(Binomial(.8, two_armed = FALSE))),
        "Binomial<single-armed>, response rate in control group: 0.8 "
    )

})
