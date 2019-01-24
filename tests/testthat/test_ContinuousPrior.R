context("ContinuousPrior                                                      ")

test_that("single point prior", {

    dist <- Normal()

    prior <- ContinuousPrior(function(x) 2*x, c(0, 1))

    expect_equal(
        bounds(prior),
        c(0, 1.))

    expect_equal(
        otsd::expectation(prior, function(x) x),
        2/3)

    n1 <- 20

    expect_equal(
        stats::integrate(
            function(z1) predictive_pdf(dist, prior, z1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)), abs.tol = .0001)$value,
        1,
        tolerance = .005)

    expect_gt(
        stats::integrate(
            function(z1) z1 * predictive_pdf(dist, prior, z1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)), abs.tol = .0001)$value,
        0)

    cprior <- condition(prior, c(.0, .5))

    expect_equal(
        bounds(cprior),
        c(.001, .5))

    expect_equal(
        stats::integrate(
            function(z1) predictive_pdf(dist, cprior, z1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)), abs.tol = .0001)$value,
        1,
        tolerance = .005)

    expect_gt(
        stats::integrate(
            function(z1) z1 * predictive_pdf(dist, prior, z1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)), abs.tol = .0001)$value,
        stats::integrate(
            function(z1) z1 * predictive_pdf(dist, cprior, z1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)), abs.tol = .0001)$value)

    delta <- 2
    z1    <- delta*sqrt(20)
    post <- posterior(dist, prior, z1, n1)

    expect_equal(
        stats::integrate(
            function(theta) post@pdf(theta),
            bounds(post)[1], bounds(post)[2], abs.tol = .0001)$value,
        1,
        tolerance = .005)

    expect_equal(
        bounds(post),
        bounds(prior))

    expect_gt(
        post@pdf(1),# good outcome should shift mass to higher values
        prior@pdf(1))

})
