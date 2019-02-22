context("GSDesign")

test_that("Group-sequential design constructor works", {
    # define an initial design
    n1      <- 25
    c1f     <-   .0
    c1e     <-  2.0
    order   <-  5L
    n2_cont <- 40
    c2_piv  <- rep(1.96, order)
    design  <- gq_design(n1, c1f, c1e, n2_cont, c2_piv, order)

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

    ess  <- integrate(ConditionalSampleSize(dist, alternative))
    cp   <- ConditionalPower(dist, alternative)
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(dist, null))
    smth <- SmoothnessN2()



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

    expect_equal(
        evaluate(smth, design),
        0.0
    )

}) # end 'group-sequential design constructor works'



test_that("Error definition works", {
    expect_error(
        GSDesign(50, 0, 2, 50, c(2, 2), c(.1, .5, .9), rep(1/3, 3))
    ) # pivots length must fit

    expect_error(
        GSDesign(50, 0, 2, 50, rep(2, 3), c(0, 1, 2), rep(1/3, 3))
    ) # x1_norm_pivots must not be scaled

    expect_error(
        GSDesign(50, 0, 2, 50, rep(2,3), c(-.5, 0, .5), c(1, 0, 1))
    ) # positive weights

}) # end 'error definition works'



test_that("GSDesign can be converted to TwoStageDesign", {
    int_pars <- GaussLegendreRule(5L)
    design1  <- GSDesign(50, 0, 2, 50, rep(2, 5), int_pars$nodes, int_pars$weights)
    design2  <- TwoStageDesign(design1)

    pow <- integrate(ConditionalPower(Normal(), PointMassPrior(.3, 1)))

    expect_equal(
        evaluate(pow, design1),
        evaluate(pow, design2)
    ) # power remains equal


    ess <- integrate(ConditionalSampleSize(Normal(), PointMassPrior(.3, 1)))

    expect_equal(
        evaluate(ess, design1),
        evaluate(ess, design2)
    ) # ess remains equal

}) # end 'GSDesign can be converted to TwoStageDesign'



test_that("Rounding works", {

    int_pars <- GaussLegendreRule(5L)
    design1  <- GSDesign(50, 0, 2, 50.2, rep(2, 5), int_pars$nodes, int_pars$weights)


    expect_equal(
        n2(design1, 1),
        50.0
    )

})
