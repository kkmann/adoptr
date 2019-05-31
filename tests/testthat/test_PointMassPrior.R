context("PointMassPrior")



test_that("single point prior", {

    dist <- Normal(two_armed = FALSE)
    prior <- PointMassPrior(.0, 1.0)

    expect_equal(
        bounds(prior),
        c(0, 0),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        adoptr::expectation(prior, function(x) x),
        0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        adoptr::expectation(prior, function(x) x + 1),
        1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    n1 <- 20

    expect_gt(
        predictive_pdf(dist, prior, 0, n1),
        predictive_pdf(dist, prior, .1, n1))

    expect_gt(
        predictive_pdf(dist, prior, 0, n1),
        predictive_pdf(dist, prior, -.1, n1))

    expect_equal(
        stats::integrate(
            function(z1) predictive_pdf(dist, prior, z1, n1),
            qnorm(.0005), qnorm(.9995), abs.tol = .0001)$value,
        1,
        tolerance = .001, scale = 1)

    expect_equal(
        stats::integrate(
            function(z1) z1 * predictive_pdf(dist, prior, z1, n1),
            qnorm(.0005), qnorm(.9995), abs.tol = .0001)$value,
        0,
        tolerance = .001, scale = 1)

    cprior <- condition(prior, c(-1 , 1))
    expect_equal(
        cprior@theta,
        .0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        cprior@mass,
        1.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    delta <- .3
    z1    <- .3*sqrt(20)
    post  <- posterior(dist, prior, z1, n1)

    expect_equal(
        post@theta,
        .0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        post@mass,
        1.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

})



test_that("multiple points prior", {

    dist  <- Normal(two_armed = FALSE)
    prior <- PointMassPrior(c(.0, .5), c(.5, .5))

    expect_equal(
        bounds(prior),
        c(0, .5),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        adoptr::expectation(prior, function(x) x),
        .25,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        adoptr::expectation(prior, function(x) x + 1),
        1.25,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    n1 <- 20
    expect_gt(
        predictive_pdf(dist, prior, 0, n1),
        predictive_pdf(dist, prior, -.1, n1))

    expect_gt(
        predictive_pdf(dist, prior, 0.5*sqrt(n1), n1),
        predictive_pdf(dist, prior, .6*sqrt(n1), n1))

    expect_equal(
        stats::integrate(
            function(z1) predictive_pdf(dist, prior, z1, n1),
            qnorm(.0005), qnorm(.9995, mean = .5*sqrt(n1)), abs.tol = .0001)$value,
        1,
        tolerance = .001, scale = 1)

    expect_gt(
        stats::integrate(
            function(z1) z1 * predictive_pdf(dist, prior, z1, n1),
            qnorm(.0005), qnorm(.9995, mean = .5*sqrt(n1)), abs.tol = .0001)$value,
        0)

    cprior <- condition(prior, c(0 , .1)) # should reduce to case above
    expect_equal(
        cprior@theta,
        .0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        cprior@mass,
        1.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    delta <- .3
    z1    <- .3*sqrt(20)
    post  <- posterior(dist, prior, z1, n1)

    expect_equal(
        post@theta,
        c(.0, .5),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_gt(
        post@mass[2],# good outcome should shift mass to higher pivot
        post@mass[1])

})



test_that("errors are defined correctly", {

    expect_error(
        PointMassPrior(theta = .3, mass = .9)) # mass must sum up to 1

    prior <- PointMassPrior(c(.5, 1.5), rep(.5,2))
    expect_error(
        condition(prior, 1)) # interval must be of length 2

    expect_error(
        condition(prior, c(0, Inf))) # interval must be finite

    expect_error(
        condition(prior, c(1, 0))) # interval[2] must be greater or equal than interval[1]

    expect_equal(
        condition(prior, c(0, 1))@mass,
        1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'errors are defined correctly'



test_that("show method returns class name", {

    prior <- PointMassPrior(.3, 1)
    expect_true(
        class(prior)[1] == capture.output(show(prior)))

}) # end 'show method returns class name'
