context("OneStageDesign")

test_that("Optimal one-stage design with point prior is computable", {
    # define an initial design
    n1      <- 60
    c1f     <-  2.0
    design  <- OneStageDesign(n1, c1f)

    # check if functions are defined correctly
    expect_equal(
        n2(design, 1.0),
        0.0
    )

    expect_equal(
        c2(design, 1.0),
        Inf
    )

    expect_equal(
        c2(design, 3.0),
        -Inf
    )

    # check if length does fit
    expect_equal(
        length(tunable_parameters(design)),
        2
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
        60.0
    )

    expect_equal(
        round(evaluate(pow, design), 3),
        0.864
    )

    expect_equal(
        round(evaluate(toer, design), 3),
        0.023
    )

    #compute optimal design

    minimize(
        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design = design,
        lower_boundary_design = update(design, c(5, -5)),
        upper_boundary_design = update(design, c(500, 5))
    ) -> d2

    d2 <- d2$design

    expect_equal(
        round(evaluate(pow, d2), 1),
        0.8
    )

    expect_equal(
        round(evaluate(toer, d2), 2),
        0.05
    )

    expect_equal(
        sign(evaluate(ess, d2) - 40),
        -1
    )

    # Check if n2 is 0 at boundaries
    expect_equal(
        n2(d2, d2@c1f),
        0.0
    )

    expect_equal(
        n2(d2, d2@c1e),
        0.0
    )

    # check if boundaries for early stopping are equal
    expect_equal(
        d2@c1f,
        d2@c1e
    )


})



test_that("error definition works", {
    expect_error(
        OneStageDesign(100)
    )

    design <- OneStageDesign(50, 2)


    expect_error(
        plot(design)
    )

}) # end 'error definition works'



test_that("OneStageDesign can be converted to TwoStageDesign", {

    design1 <- OneStageDesign(87.21, 1.96)
    design2 <- TwoStageDesign(design1, order = 7)

    pow <- expected(ConditionalPower(Normal(two_armed = FALSE), PointMassPrior(.3, 1)))

    expect_true(
        length(design2@x1_norm_pivots) == 7
    )

    expect_true(
        length(design2@n2_pivots) == 7
    )

    expect_equal(
        evaluate(pow, design1),
        evaluate(pow, design2),
        tolerance = .001
    )

    expect_equal(
        evaluate(pow, design1),
        .8,
        tolerance = .01
    ) # power works for OneStageDesign

}) # end 'OneStageDesign can be converted to TwoStageDesign'



test_that("show method returns design name", {
    design  <-  OneStageDesign(90, 2.0)

    expect_equal(
        cat(class(design)[1]),
        show(design)
    )
}) # end 'show method returns design name'

