context("test sample size")

design <- TwoStageDesign(25, 0, 2, rep(40.0, 5), rep( 2, 5), 5L)
dist   <- Normal()
z1     <- seq(-1, 3, .1)


test_that("conditional sample size maps to actual sample size", {

    prior  <- ContinuousPrior(function(x) rep(1/10, length(x)), c(-4, 6))

    css    <- ConditionalSampleSize()

    expect_equal(
        evaluate(css, design, z1),
        n(design, z1)
    )


}) # end 'conditional sample size maps to actual sample size'



context("Test Conditional Power")

test_that("Conditional Power is monotonous", {

    # Define two simple designs
    design1 <- TwoStageDesign(25, 0, 2, rep(40.0, 5), seq(2.0, 0.0, length.out = 5))
    design2 <- TwoStageDesign(25, 0, 2, rep(60.0, 5), seq(2.0, 0.0, length.out = 5))

    dist        <- Normal()
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)

    cp <- ConditionalPower(dist, alternative)

    # Conditional Power is monotonously increasing in n
    expect_gt(
        evaluate(cp, design2, 1),
        evaluate(cp, design1, 1)
    )

    # Conditional Power is monotonously increasing in z1
    expect_gt(
        evaluate(cp, design1, 1.5),
        evaluate(cp, design1, 0.5)
    )

}) # end "Conditional Power is monotonous"



test_that("Conditional power has correct values outside continuation region",{
    design1 <- TwoStageDesign(25, 0, 2, rep(40.0, 5), seq(2.0, 0.0, length.out = 5))
    cp      <- ConditionalPower(Normal(), PointMassPrior(.4, 1))

    expect_equal(
        evaluate(cp, design1, -1),
        0.0
    )

    expect_equal(
        evaluate(cp, design1, 3),
        1.0
    )
}) # end 'Conditional power has correct values outside continuation region'


context("Test class ConditionalScore")

test_that("Conditional scores are vectorized in z1", {

    css    <- ConditionalSampleSize()
    cp     <- ConditionalPower(Normal(), PointMassPrior(.4, 1))

    expect_length(
        evaluate(css, design, z1),
        length(z1)
    )

    expect_length(
        evaluate(cp, design, z1),
        length(z1)
    )

}) # end 'Conditional scores are vectorized in z1'
