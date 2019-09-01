context("GroupSequentialDesign")



test_that("Group-sequential design constructor works", {

    order   <-  5L
    design  <-  GroupSequentialDesign(25, 0.0, 2.0, 40.0, 1.96, order)

    # check if functions are defined correctly
    expect_equal(
        n2(design, 1.0),
        40.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        c2(design, 1.0),
        1.96,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    # check if length does fit
    expect_true(
        length(tunable_parameters(design)) == order + 4)

    # check if key figures can be computed
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    dist        <- Normal(two_armed = FALSE)
    ess         <- ExpectedSampleSize(dist, alternative)
    pow         <- Power(dist, alternative)
    toer        <- Power(dist, null)

    expect_equal(
        evaluate(ess, design),
        44.1,
        tolerance = 1e-1, scale = 1)

    expect_equal(
        evaluate(pow, design),
        0.842,
        tolerance = 1e-3, scale = 1)

    expect_equal(
        evaluate(toer, design),
        0.035,
        tolerance = 1e-3, scale = 1)

})  # end 'group-sequential design constructor works'





test_that("GSDesign can be converted to TwoStageDesign", {

    design1 <- GroupSequentialDesign(50, 0, 2, 50, rep(2, 5))
    design2 <- TwoStageDesign(design1)
    pow     <- Power(Normal(), PointMassPrior(.3, 1))
    ess     <- ExpectedSampleSize(Normal(), PointMassPrior(.3, 1))

    expect_equal(
        evaluate(pow, design1),
        evaluate(pow, design2),
        tolerance = 1e-3, scale = 1)

    expect_equal(
        evaluate(ess, design1),
        evaluate(ess, design2),
        tolerance = 1e-3, scale = 1)

}) # end 'GSDesign can be converted to TwoStageDesign'



test_that("Rounding works", {

    expect_equal(
        n2(GroupSequentialDesign(50, 0, 2, 50.2, rep(2, 5)), 1),
        50.0,
        tolerance = 1e-6, scale = 1
    )

})



test_that("show method", {

    expect_equal(
        paste0(capture.output(show(GroupSequentialDesign(25, 0.0, 2.0, 40.0, 1.96, 5L))), collapse = "\n\r"),
        "GroupSequentialDesign<\n\r\r      x1    c2   n\n\r\r    -0.00   Inf   25\n\r\r     0.09  1.96   65\n\r\r     0.46  1.96   65\n\r\r     1.00  1.96   65\n\r\r     1.54  1.96   65\n\r\r     1.91  1.96   65\n\r\r     2.00  -Inf   25>\n\r\r "
    )

})
