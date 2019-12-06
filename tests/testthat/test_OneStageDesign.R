context("OneStageDesign")



test_that("Optimal one-stage design with point prior is computable", {

    # define an initial design
    n1      <- 60
    c1f     <-  2.0
    design  <- OneStageDesign(n1, c1f)

    # check if functions are defined correctly
    expect_equal(
        n2(design, 1.0),
        0.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_true(
        c2(design, 1.0) == Inf)

    expect_true(
        c2(design, 3.0) == -Inf)

    # check if length does fit
    expect_true(
        length(tunable_parameters(design)) == 2)

    # check if key figures can be computed
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    dist        <- Normal(two_armed = FALSE)
    ess         <- ExpectedSampleSize(dist, alternative)
    pow         <- Power(dist, alternative)
    toer        <- Power(dist, null)

    expect_equal(
        evaluate(ess, design),
        60.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(pow, design),
        0.864,
        tolerance = 1e-3, scale = 1)

    expect_equal(
        evaluate(toer, design),
        0.023,
        tolerance = 1e-3, scale = 1)

    #compute optimal design

    d2 <- minimize(
        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        design
    )

    d2 <- d2$design

    expect_equal(
        round(evaluate(pow, d2), 1),
        0.8,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        round(evaluate(toer, d2), 2),
        0.05,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_true(
        sign(evaluate(ess, d2) - 40) == -1)

    # Check if n2 is 0 at boundaries
    expect_equal(
        n2(d2, d2@c1f),
        0.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        n2(d2, d2@c1e),
        0.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    # check if boundaries for early stopping are equal
    expect_equal(
        d2@c1f,
        d2@c1e,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

})



test_that("error definition works", {

    expect_error(
        OneStageDesign(100))

    design <- OneStageDesign(50, 2)

    expect_error(
        plot(design))

}) # end 'error definition works'



test_that("OneStageDesign can be converted to TwoStageDesign", {

    design1 <- OneStageDesign(87.21, 1.96)
    design2 <- TwoStageDesign(design1, order = 7)
    pow     <- Power(Normal(two_armed = FALSE), PointMassPrior(.3, 1))

    expect_true(length(design2@x1_norm_pivots) == 7)
    expect_true(length(design2@n2_pivots) == 7)

    expect_equal(
        evaluate(pow, design1),
        evaluate(pow, design2),
        tolerance = 1e-6, scale = 1)

    expect_equal(
        evaluate(pow, design1, optimization = TRUE),
        .8,
        tolerance = 1e-5, scale = 1)

    expect_lte(
        evaluate(pow, design1, optimization = FALSE),
        .8)


}) # end 'OneStageDesign can be converted to TwoStageDesign'



test_that("show method returns design name", {

    expect_equal(
        capture.output(show(OneStageDesign(90, 2.0))),
        "OneStageDesign<n=90;c=2.00> "
    )

})
