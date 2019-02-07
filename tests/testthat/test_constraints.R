context("constraint specifications                                            ")

# TODO: values are ad-hoc, come up with a way to actually verify them!

test_that("UnconditionalConstraints", {

    # create dummy design
    design <- gq_design(25, 0, 2, rep(40.0, 5), rep(1.96, 5), 5L)

    # create power as IntegralScore
    pow <- integrate(ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1)))

    # construct actual constraint
    cnstr <- pow >= 0.8

    # see if it evaluates to the right value
    expect_equal(
        evaluate(cnstr, design), -0.0415, tolerance = .001)

    # check other direction
    toer <- integrate(ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.0, 1)))
    expect_equal(
        evaluate(toer <= .05, design), -.0153, tolerance = .001)

})

test_that("ConditionalConstraints", {

    # create dummy design
    design <- gq_design(25, 0, 2, rep(40.0, 5), rep(1.96, 5), 5L)

    # create power as IntegralScore
    cp <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1))

    # construct actual constraint
    cnstr <- cp >= 0.8

    # see if it evaluates to the right value
    expect_equal(
        evaluate(cnstr, design, .8), 0.0844, tolerance = .001)

    # check other direction
    ctoer <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.0, 1))
    expect_equal(
        evaluate(ctoer <= .05, design, .8), -.0250, tolerance = .001)

})


test_that("ConditionalConstraints", {

    # create dummy design
    design <- gq_design(25, 0, 2, rep(40.5, 5), rep(1.96, 5), 5L)

    # create power as IntegralScore
    cp <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1))

    # construct a constraint set and see if it is at least of the right length
    expect_equal(
        length(
            evaluate(subject_to(cp >= .6, cp >= .5), design)
        ),
        10)

    # Compute conditional sample size
    css <- ConditionalSampleSize(Normal(), PointMassPrior(.4, 1))

    # Use non-rounded values
    expect_equal(
        evaluate(css, design, 1),
        65.5
    )

    # Use rounded values
    design@rounded <- TRUE
    expect_equal(
        evaluate(css, design, 1),
        65.0
    )

})



test_that("subject_to throws correct error", {

    expect_error(subject_to(1))

})
