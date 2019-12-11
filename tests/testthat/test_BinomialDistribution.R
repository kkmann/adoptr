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


# Let us check if our tests are working correctly
r_c    <<- .2
r_e    <<- .3
r_diff <<- r_e - r_c
dist   <<- Binomial(r_c)
alpha  <<- 0.025
beta   <<- 0.2
n      <<- ceiling(rpact::getSampleSizeRates(alpha = alpha, beta = beta,
                                             pi1 = r_e, pi2 = r_c)$nFixed1)
c      <<- qnorm(1 - alpha)



test_that("Approximation works in cdf", {
    iters   <- 1e4
    pval_H0 <- rep(0, iters)
    pval_H1 <- rep(0, iters)

    for(i in 1:iters) {
        x_c   <- stats::rbinom(1, n, r_c)
        x_e   <- stats::rbinom(1, n, r_e)
        p_c   <- x_c / n
        p_e   <- x_e / n
        p_0   <- (x_c + x_e) / (2 * n)
        test_statistic_H1 <- sqrt(n / 2) * (p_e - p_c) / sqrt(p_0 * (1 - p_0))
        pval_H1[i] <- 1 - cumulative_distribution_function(dist, test_statistic_H1, n, 0)

        x_e_2 <- stats::rbinom(1, n, r_c)
        p_e_2 <- x_e_2 / n
        p_0_2 <- (x_c + x_e_2) / (2 * n)
        test_statistic_H0 <- sqrt(n / 2) * (p_e_2 - p_c) / sqrt(p_0_2 * (1 - p_0_2))
        pval_H0[i] <- 1 - cumulative_distribution_function(dist, test_statistic_H0, n, 0)
    }

    expect_equal(
        mean(pval_H1 <= alpha),
        1 - beta,
        tolerance = 1 / sqrt(iters), scale = 1
    )

    expect_equal(
        mean(pval_H0 <= alpha),
        alpha,
        tolerance = 1 / sqrt(iters), scale = 1
    )

})


test_that("Approximation works in cdf", {
    design <- OneStageDesign(n, c)
    H_0    <- PointMassPrior(0, 1)
    H_1    <- PointMassPrior(r_diff, 1)

    expect_equal(
        evaluate(Power(dist, H_1), design),
        1 - beta,
        tolerance = 1 / sqrt(iters), scale = 1
    )

    expect_equal(
        evaluate(Power(dist, H_0), design),
        alpha,
        tolerance = 1 / sqrt(iters), scale = 1
    )

})




test_that("Approximation works in simulate", {
    iters  <- 1e6
    test_statistic_H0 <- simulate(object = dist, nsim = iters, n = n, theta = 0, seed = 42)
    test_statistic_H1 <- simulate(object = dist, nsim = iters, n = n, theta = r_diff, seed = 42)
    mean(test_statistic_H1 > c)

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
