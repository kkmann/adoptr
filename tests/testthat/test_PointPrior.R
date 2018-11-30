context("PointMassPrior                                                       ")

test_that("single point priors", {

    prior <- PointMassPrior(.0, 1.0)

    expect_equal(
        bounds(prior),
        c(0, 0))

    expect_equal(
        otsd::expectation(prior, function(x) x),
        0)

    expect_equal(
        otsd::expectation(prior, function(x) x + 1),
        1)

    n1    <- 20

    expect_gt(
        predictive_pdf(prior, 0, n1),
        predictive_pdf(prior, .1, n1))

    expect_gt(
        predictive_pdf(prior, 0, n1),
        predictive_pdf(prior, -.1, n1))

    expect_equal(
        stats::integrate(
            function(z1) predictive_pdf(prior, z1, n1),
            qnorm(.0005), qnorm(.9995), abs.tol = .0001)$value,
        1,
        tolerance = .005)

    expect_equal(
        stats::integrate(
            function(z1) z1 * predictive_pdf(prior, z1, n1),
            qnorm(.0005), qnorm(.9995), abs.tol = .0001)$value,
        0,
        tolerance = .005)

    cprior <- condition(prior, c(-1 , 1))

    expect_equal(
        cprior@theta,
        .0)

    expect_equal(
        cprior@mass,
        1.0)

    delta <- .3
    z1    <- .3*sqrt(20)
    post <- posterior(prior, z1, n1)

    expect_equal(
        post@theta,
        .0)

    expect_equal(
        post@mass,
        1.0)

})
