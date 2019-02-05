context("GSDesign")

test_that("Optimal group-sequential design with point prior is computable", {
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
    smth <- integrate(SmoothnessN2(dist))



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

    #compute optimal design

    minimize(
        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design = design,
        lower_boundary_design = update(design, c(10, -1, 1, 2, numeric(order) - 5)),
        upper_boundary_design = update(design, c(50, 1, 4, 50, numeric(order) + 5))
    ) ->
        d2


    expect_equal(
        round(evaluate(pow, d2), 1),
        0.8
    )

    expect_equal(
        round(evaluate(toer, d2), 2),
        0.05
    )

    expect_equal(
        sign(evaluate(ess, d2) - 32),
        -1
    )

    # Check if n2 is equal at boundaries
    expect_equal(
        n2(d2, d2@c1f),
        n2(d2, d2@c1e)
    )


})
