context("Test ConditionalSampleSize                                          ")

test_that("conditional sample size maps to actual sample size", {

    design <- GQDesign(25, 0, 2, rep(40.0, 5), rep( 1.96, 5), 5)

    dist <- Normal()
    z1   <- seq(-1, 3, .1)

    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.6, 1)

    css0 <- ConditionalSampleSize(dist, null)
    css1 <- ConditionalSampleSize(dist, alternative)

    expect_equal(
        evaluate(css0, design, z1),
        n(design, z1))

    expect_equal(
        evaluate(css0, design, z1),
        evaluate(css1, design, z1))

    # basic checks on integral score
    ess0 <- integrate(css0)
    ess1 <- integrate(css1)

    expect_gt(
        evaluate(ess0, design),
        evaluate(ess1, design)
    )

    expect_gt(
        n(design, 1),
        evaluate(ess1, design)
    )

    expect_gt(
        evaluate(ess1, design),
        n1(design)
    )

})



test_that("conditional score arithmetic works", {

    design <- GQDesign(25, 0, 2, rep(40.0, 5), rep( 1.96, 5), 5)

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
