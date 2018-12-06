context("constraint specifications                                            ")

# TODO: values are ad-hoc, come up with a way to actually verify them!

test_that("UnconditionalConstraints", {

    # create dummy design
    design <- GQDesign(25, 0, 2, rep(40.0, 5), rep( 1.96, 5), 5)

    # create power as IntegralScore
    pow <- integrate(ConditionalPower(Normal(), PointMassPrior(.4, 1)))

    # construct actual constraint
    cnstr <- pow >= 0.8

    # see if it evaluates to the right value
    expect_equal(
        evaluate(cnstr, design), -0.0415, tolerance = .001)

    # check other direction
    toer <- integrate(ConditionalPower(Normal(), PointMassPrior(.0, 1)))
    expect_equal(
        evaluate(toer <= .05, design), -.0153, tolerance = .001)

})

test_that("ConditionalConstraints", {

    # create dummy design
    design <- GQDesign(25, 0, 2, rep(40.0, 5), rep( 1.96, 5), 5)

    # create power as IntegralScore
    cp <- ConditionalPower(Normal(), PointMassPrior(.4, 1))

    # construct actual constraint
    cnstr <- cp >= 0.8

    # see if it evaluates to the right value
    expect_equal(
        evaluate(cnstr, design, .8), 0.0844, tolerance = .001)

    # check other direction
    ctoer <- ConditionalPower(Normal(), PointMassPrior(.0, 1))
    expect_equal(
        evaluate(ctoer <= .05, design, .8), -.0250, tolerance = .001)

})
