context("GroupSequentialDesign")

test_that("Group-sequential design constructor works", {
    # define an initial design
    order   <-  5L
    design  <-  GroupSequentialDesign(25, 0.0, 2.0, 40.0, 1.96, order)

    # check if functions are defined correctly
    expect_equal(
        n2(design, 1.0),
        40.0
    )

    expect_equal(
        c2(design, 1.0),
        1.96
    )

    # check if length does fit
    expect_equal(
        length(tunable_parameters(design)),
        order + 4
    )

    # check if key figures can be computed
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)

    dist <- Normal(two_armed = FALSE)

    ess  <- expected(ConditionalSampleSize(dist, alternative))
    cp   <- ConditionalPower(dist, alternative)
    pow  <- expected(cp)
    toer <- expected(ConditionalPower(dist, null))


    expect_equal(
        round(evaluate(ess, design), 1),
        44.1
    )

    expect_equal(
        round(evaluate(pow, design), 3),
        0.842
    )

    expect_equal(
        round(evaluate(toer, design), 3),
        0.035
    )

})  # end 'group-sequential design constructor works'





test_that("GSDesign can be converted to TwoStageDesign", {

    design1  <- GroupSequentialDesign(50, 0, 2, 50, rep(2, 5))
    design2  <- TwoStageDesign(design1)

    pow <- expected(ConditionalPower(Normal(), PointMassPrior(.3, 1)))

    expect_equal(
        evaluate(pow, design1),
        evaluate(pow, design2)
    ) # power remains equal


    ess <- expected(ConditionalSampleSize(Normal(), PointMassPrior(.3, 1)))

    expect_equal(
        evaluate(ess, design1),
        evaluate(ess, design2)
    ) # ess remains equal

}) # end 'GSDesign can be converted to TwoStageDesign'



test_that("Rounding works", {

    design1  <- GroupSequentialDesign(50, 0, 2, 50.2, rep(2, 5))


    expect_equal(
        n2(design1, 1),
        50.0
    )

})


test_that("show method returns design name", {
    design  <-  GroupSequentialDesign(25, 0.0, 2.0, 40.0, 1.96, 5L)

    expect_equal(
        cat(class(design)[1]),
        show(design)
    )
}) # end 'show method returns design name'

