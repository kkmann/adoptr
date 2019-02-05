context("test sample size")

test_that("conditional sample size maps to actual sample size", {

    design <- gq_design(25, 0, 2, rep(40.0, 5), rep( 2, 5), 5L)
    dist   <- Normal()
    z1     <- seq(-1, 3, .1)
    prior  <- ContinuousPrior(function(x) rep(1/10, length(x)), c(-4, 6))

    css    <- ConditionalSampleSize(dist, prior)

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


context("Test Conditional Power")

test_that("Conditional Power is monotonous", {
    # Define two simple designs
    design1 <- gq_design(25, 0, 2, rep(40.0, 5), seq(2.0, 0.0, length.out = 5), 5L)
    design2 <- gq_design(25, 0, 2, rep(60.0, 5), seq(2.0, 0.0, length.out = 5), 5L)

    dist <- Normal()
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)

    cp <- ConditionalPower(dist, alternative)

    # Conditional Power is monotonously increasing in n
    expect_gt(evaluate(cp, design2, 1),
              evaluate(cp, design1, 1))

    # Conditional Power is monotonously increasing in z1
    expect_gt(evaluate(cp, design1, 1.5),
              evaluate(cp, design1, 0.5))

}) # end "Conditional Power is monotonous"



test_that("Conditional power has correct values outside continuation region",{
    design1 <- gq_design(25, 0, 2, rep(40.0, 5), seq(2.0, 0.0, length.out = 5), 5L)
    cp      <- ConditionalPower(Normal(), PointMassPrior(.4, 1))

    expect_equal(
        evaluate(cp, design1, -1),
        0)

    expect_equal(
        evaluate(cp, design1, 3),
        1)
}) # end 'Conditional power has correct values outside continuation region'


context("Test vectorization")

test_that("Conditional scores are vectorized in z1", {
    design <- gq_design(25, 0, 2, rep(40.0, 5), rep( 1.96, 5), 5L)
    z1     <- seq(-1, 3, .1)
    css    <- ConditionalSampleSize(Normal(), PointMassPrior(.4, 1))
    cp     <- ConditionalPower(Normal(), PointMassPrior(.4, 1))

    expect_length(
        evaluate(css, design, z1),
        length(z1))

    expect_length(
        evaluate(cp, design, z1),
        length(z1))

}) # end 'Conditional scores are vectorized in z1'
