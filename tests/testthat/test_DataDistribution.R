context("Normal data distribution")

test_that("Constructor works", {

    expect_true(
        Normal()@two_armed
    )

    expect_true(
        !Normal(two_armed = FALSE)@two_armed
    )

}) # end 'constructor works'



test_that("pdf corresponds to dnorm", {

    x     <- seq(-3, 3, by = .1)
    n     <- 10
    theta <- .5
    expect_equal(
        probability_density_function(
            Normal(two_armed = FALSE),
            x, n, theta
        ),
        dnorm(x, sqrt(n) * theta, 1)
    )

})



test_that("cdf corresponds to pnorm", {

    x     <- seq(-3, 3, by = .1)
    n     <- 10
    theta <- .5
    expect_equal(
        cumulative_distribution_function(
            Normal(),
            x, n, theta
        ),
        pnorm(x, sqrt(n) * theta / sqrt(2), 1)
    )

})



test_that("quantile corresponds to qnorm", {

    probs <- seq(0.1, 1, by = .1)
    n     <- 50
    theta <- .25
    expect_equal(
        quantile(
            Normal(),
            probs, n, theta
        ),
        qnorm(probs, sqrt(n) * theta / sqrt(2), 1)
    )

    expect_warning(
        quantile(Normal(), -1, n, theta)
    )

    expect_equal(
        suppressWarnings(quantile(Normal(), -1, n, theta)),
        NaN
    )

})



test_that("simulate respects seed", {

    expect_equal(
        simulate(Normal(), 10, 25, .5, seed = 42),
        simulate(Normal(), 10, 25, .5, seed = 42)
    )

    set.seed(42)

    expect_true(
        all(simulate(Normal(), 10, 25, -.5) != simulate(Normal(), 10, 25, -.5))
    )
}) # end 'simulate respects seed'



test_that("show method returns class name", {
    dist <- Normal()

    expect_equal(
        cat(class(dist)[1]),
        show(dist)
    )
}) # end 'show method returns class name'

