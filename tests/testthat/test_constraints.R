context("constraint specifications                                            ")




design <- TwoStageDesign(25, 0, 2, 40.5, 1.96, 5L)



test_that("UnconditionalConstraints", {

    pow   <- Power(Normal(two_armed = FALSE), PointMassPrior(.4, 1))
    toer  <- Power(Normal(two_armed = FALSE), PointMassPrior(.0, 1))
    cnstr <- pow >= 0.8
    pow_true <-  mean(adoptr::simulate
                      (design, nsim = 10^6, dist = Normal(two_armed = FALSE),
                          theta = .4, seed = 42)$reject)

    # see if it evaluates to the right value
    expect_equal(
        evaluate(cnstr, design),
        (.8 - pow_true),
        tolerance = 1e-3, scale = 1)

    # compute true value
    toer_true <-  mean(adoptr::simulate
                      (design, nsim = 10^6, dist = Normal(two_armed = FALSE),
                          theta = .0, seed = 142)$reject)

    expect_equal(
        evaluate(toer <= .05, design),
        (toer_true - .05),
        tolerance = 1e-3, scale = 1)

    # Check syntax
    expect_true(
        evaluate(subject_to(pow >= .8), design) == evaluate(subject_to(.8 <= pow), design))

    expect_true(
        evaluate(subject_to(.05 >= toer), design) == evaluate(subject_to(toer <= .05), design))

})



test_that("ConditionalConstraints", {

    # create conditional power
    cp <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1))
    cnstr <- cp >= 0.8

    expect_equal(
        evaluate(cnstr, design, .8),
        0.0844,
        tolerance = 1e-3, scale = 1)

    # check other direction
    ctoer <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.0, 1))

    expect_equal(
        evaluate(ctoer <= .05, design, .8),
        -.0250,
        tolerance = 1e-3, scale = 1)

})



test_that("ConditionalConstraints", {

    cp <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1))
    css <- ConditionalSampleSize()

    # construct a constraint set and see if it is at least of the right length
    expect_true(
        length(evaluate(subject_to(cp >= .6, cp >= .5), design)) == 10)

    # Use non-rounded values
    expect_equal(
        evaluate(css, design, 1, optimization = TRUE),
        65.5,
        tolerance = 1e-3, scale = 1)

    # Use rounded values
    expect_equal(
        evaluate(css, design, 1),
        65.0,
        tolerance = 1e-3, scale = 1)

    # Check syntax
    expect_true(all(
        evaluate(subject_to(cp >= .8), design) == evaluate(subject_to(.8 <= cp), design)))

    expect_true(all(
        evaluate(subject_to(css <= 500), design) == evaluate(subject_to(500 >= css), design)))

})



test_that("subject_to throws correct error", {

    expect_error(subject_to(1))

}) # end 'subject_to throws correct error'



test_that("score vs score inequalities", {

    cp    <- ConditionalPower(Normal(), PointMassPrior(.28, 1))
    ctoer <- ConditionalPower(Normal(), PointMassPrior(.0, 1))
    pow   <- expected(cp, Normal(), PointMassPrior(.28, 1))
    toer  <- expected(ctoer, Normal(), PointMassPrior(.0, 1))

    expect_true(all(
        evaluate(subject_to(ctoer <= cp), design) == evaluate(subject_to(cp >= ctoer), design)))

    expect_true(
        evaluate(subject_to(toer <= pow), design) == evaluate(subject_to(pow >= toer), design))

}) # end 'score vs score inequalities'



test_that("show method returns class name", {

    cp <- ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.4, 1))
    cnstr <- cp >= 0.8

    expect_true(
        capture.output(show(cnstr)) == class(cnstr)[1])

})
