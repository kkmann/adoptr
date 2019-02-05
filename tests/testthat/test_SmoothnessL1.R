context("Smoothness")

test_that("SmoothnessL1 defines a norm", {
    order <- 5L

    design0 <- gq_design(25, 0, 2, rep(0.0, order), rep(2, order), order)


    design1 <- gq_design(25, 0, 2, rep(40.0, order), rep(2, order), order)

    design2 <- design1
    design2@n2_pivots <- 2 * design1@n2_pivots

    design3 <- design0
    design3@n2_pivots <- design1@n2_pivots + design2@n2_pivots

    smth <- SmoothnessL1()

    expect_equal(
        evaluate(smth, design0),
        0.0
    )

    expect_equal(
        evaluate(smth, design2),
        2 * evaluate(smth, design1)
    )

    expect_equal( # everything is positive!
        evaluate(smth, design3),
        evaluate(smth, design1) + evaluate(smth, design2)
    )
}) # end 'SmoothnessL1 defines a norm'

test_that("Specific implementation of SmoothnessL1 yields similar results", {
    design <- gq_design(25, 0, 2, seq(60.0, 10.0, length.out = order), rep(2, order), order)
    smth   <- SmoothnessL1()

    expect_equal(
        evaluate(smth, design),
        evaluate(smth, design, specific = FALSE)
    )

    # Define something weird
    design2 <- gq_design(25, 0, 2, c(20, 72, 15, 3, 80), rep(2, order), order)

    expect_equal(
        evaluate(smth, design2),
        evaluate(smth, design2, specific = FALSE),
        tolerance = 1.0
    )


}) # end 'specific implementation of SmoothnessL1 yields similar results'
