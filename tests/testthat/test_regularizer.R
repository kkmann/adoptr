context("Regularization")

test_that("AverageN2 defines a norm", {
    order <- 5L

    design0 <- gq_design(25, 0, 2, rep(0.0, order), rep(2, order), order)


    design1 <- gq_design(25, 0, 2, rep(40.0, order), rep(2, order), order)

    design2 <- design1
    design2@n2_pivots <- 2 * design1@n2_pivots

    design3 <- design0
    design3@n2_pivots <- design1@n2_pivots + design2@n2_pivots

    avn2 <- AverageN2()

    expect_equal(
        evaluate(avn2, design0),
        0.0
    )

    expect_equal(
        evaluate(avn2, design2),
        2 * evaluate(avn2, design1)
    )

    expect_equal( # everything is positive!
        evaluate(avn2, design3),
        evaluate(avn2, design1) + evaluate(avn2, design2)
    )
}) # end 'AverageN2 defines a norm'



test_that("Specific implementation of AverageN2 yields similar results", {
    order = 9L
    design <- gq_design(25, 0, 2, seq(60.0, 10.0, length.out = order), rep(2, order), order)
    avn2   <- AverageN2()

    expect_equal(
        evaluate(avn2, design),
        evaluate(avn2, design, specific = FALSE)
    )

    # Define something weird
    order = 5L
    design2 <- gq_design(25, 0, 2, c(20, 72, 15, 3, 80), rep(2, order), order)

    expect_equal(
        evaluate(avn2, design2),
        evaluate(avn2, design2, specific = FALSE),
        tolerance = 1.0
    )


}) # end 'specific implementation of AverageN2 yields similar results'



test_that("SmoothnessN2 is defined sensitive", {
    order = 7L
    design <- gq_design(25, 0, 2, rep(40.0, order), rep(2, order), order)
    smth   <- SmoothnessN2()

    expect_equal(
        evaluate(smth, design),
        0.0
    )

    expect_equal(
        evaluate(smth, design, specific = FALSE),
        0.0
    )


    design1 <- gq_design(25, 0, 2, seq(60.0, 10.0, length.out = order), rep(2, order), order)
    design2 <- design1
    design2@n2_pivots[1] <- 70.0

    expect_gt(
        evaluate(smth, design2),
        evaluate(smth, design1)
    )

    expect_gt(
        evaluate(smth, design2, specific = FALSE),
        evaluate(smth, design1, specific = FALSE)
    )


}) # end 'SmoothnessN2 is defined sensitive'



test_that("PenaltyN1 works", {
    order = 7L
    design1 <- gq_design(25, 0, 2, rep(40.0, order), rep(2, order), order)
    design2 <- gq_design(55, 0, 2, rep(40.0, order), rep(2, order), order)

    pN1 <- PenaltyN1()

    expect_gt(
        evaluate(pN1, design2),
        evaluate(pN1, design1)
    )

}) # end 'PenaltyN1 works'
