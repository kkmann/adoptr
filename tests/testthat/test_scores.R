context("test sample size")

test_that("conditional sample size maps to actual sample size", {

    design <- gq_design(25, 0, 2, rep(40.0, 5), rep( 2, 5), 5L)
    dist   <- Normal()
    z1     <- seq(-1, 3, .1)
    prior  <- ContinuousPrior(function(x) rep(1/10, length(x)), c(-4, 6))

    css    <<- ConditionalSampleSize(dist, prior)

    expect_equal(
        evaluate(css, design, z1),
        n(design, z1)
    )

}) # end 'conditional sample size maps to actual sample size'




test_that("conditional score arithmetic works", {

    design <- gq_design(25, 0, 2, rep(40.0, 5), rep( 1.96, 5), 5L)

    dist <- Normal()
    z1   <- seq(-1, 3, .1)

    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.6, 1)

    css0 <- ConditionalSampleSize(dist, null)
    css1 <- ConditionalSampleSize(dist, alternative)

    expect_equal(
        evaluate(2*css0 + css1, design, z1),
        2*evaluate(css0, design, z1) + evaluate(css1, design, z1))

    # more test cases!

})
