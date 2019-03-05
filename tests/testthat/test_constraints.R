context("constraint specifications                                            ")

# create dummy design
design <- TwoStageDesign(25, 0, 2, 40.5, 1.96, 5L)


test_that("UnconditionalConstraints", {

    # create power as IntegralScore
    pow <- expected(ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1)))

    # construct actual constraint
    cnstr <- pow >= 0.8

    # compute true value
    pow_true <-  mean(adoptr::simulate
                      (design, nsim = 10^6, dist = Normal(two_armed = FALSE),
                          theta = .4, seed = 42)$reject)

    # see if it evaluates to the right value
    expect_equal(
        evaluate(cnstr, design), (.8 - pow_true), tolerance = .001)

    # check other direction
    toer <- expected(ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.0, 1)))
    # compute true value
    toer_true <-  mean(adoptr::simulate
                      (design, nsim = 10^6, dist = Normal(two_armed = FALSE),
                          theta = .0, seed = 142)$reject)

    expect_equal(
        evaluate(toer <= .05, design), (toer_true - .05), tolerance = .001)


    # Check syntax
    expect_equal(
        evaluate(subject_to(pow >= .8), design),
        evaluate(subject_to(.8 <= pow), design)
    )


    # Check syntax
    expect_equal(
        evaluate(subject_to(.05 >= toer), design),
        evaluate(subject_to(toer <= .05), design)
    )


})

test_that("ConditionalConstraints", {


    # create conditional power
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


    # create conditional power
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
        evaluate(css, design, 1, optimization = TRUE),
        65.5
    )

    # Use rounded values
    expect_equal(
        evaluate(css, design, 1),
        65.0
    )


    # Check syntax
    expect_equal(
        evaluate(subject_to(cp >= .8), design),
        evaluate(subject_to(.8 <= cp), design)
    )


    # Check syntax
    expect_equal(
        evaluate(subject_to(css <= 500), design),
        evaluate(subject_to(500 >= css), design)
    )

})



test_that("subject_to throws correct error", {

    expect_error(subject_to(1))

}) # end 'subject_to throws correct error'



test_that("score vs score inequalities", {

    # create conditional scores
    cp    <- ConditionalPower(Normal(), PointMassPrior(.28, 1))
    ctoer <- ConditionalPower(Normal(), PointMassPrior(.0, 1))

    expect_equal(
        evaluate(subject_to(ctoer <= cp), design),
        evaluate(subject_to(cp >= ctoer), design)
    )


    # create unconditional scores
    pow  <- expected(cp)
    toer <- expected(ctoer)

    expect_equal(
        evaluate(subject_to(toer <= pow), design),
        evaluate(subject_to(pow >= toer), design)
    )

}) # end 'score vs score inequalities'



test_that("show method returns class name", {
    cp <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1))
    cnstr <- cp >= 0.8

    expect_equal(
        show(cnstr),
        cat(class(cnstr)[1])
    )
}) # end 'show method returns class name'

