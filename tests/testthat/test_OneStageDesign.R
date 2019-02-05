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

    ess  <- integrate(ConditionalSampleSize(dist, alternative))
    cp   <- ConditionalPower(dist, alternative)
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(dist, null))
    smth <- integrate(SmoothnessN2(dist))



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

    expect_equal(
        evaluate(smth, design),
        0.0
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
