context("Regularization")



test_that("AverageN2 defines a norm", {

    order <- 5L
    design0 <- TwoStageDesign(25, 0, 2, 0.0, 2.0, order)
    design1 <- TwoStageDesign(25, 0, 2, 40.0, 2.0, order)

    design2 <- design1
    design2@n2_pivots <- 2 * design1@n2_pivots

    design3 <- design0
    design3@n2_pivots <- design1@n2_pivots + design2@n2_pivots

    avn2 <- AverageN2()

    expect_equal(
        evaluate(avn2, design0),
        0.0,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(avn2, design2),
        2 * evaluate(avn2, design1),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal( # everything is positive!
        evaluate(avn2, design3),
        evaluate(avn2, design1) + evaluate(avn2, design2),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'AverageN2 defines a norm'



test_that("Specific implementation of AverageN2 yields similar results", {

    order  <- 9L # integrates polynomials up to order 2*9 - 1 exactly,
                 # linear must be exact!
    design <- TwoStageDesign(25, 0, 2, seq(60.0, 10.0, length.out = order), rep(2, order))
    avn2   <- AverageN2()
    expect_equal(
        evaluate(avn2, design),
        evaluate(avn2, design, optimization = FALSE),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    # define something weird - accept difference
    order = 5L
    design2 <- TwoStageDesign(25, 0, 2, c(20, 72, 15, 3, 80), rep(2, order))
    expect_equal(
        evaluate(avn2, design2),
        evaluate(avn2, design2, optimization = TRUE),
        tolerance = .5, scale = 1)

}) # end 'specific implementation of AverageN2 yields similar results'



test_that("N1 works", {

    order   <- 7L
    design1 <- TwoStageDesign(25.75, 0, 2, 40.0, 2.0, order)
    pN1 <- N1()

    expect_equal(
        evaluate(pN1, design1),
        26,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(pN1, design1, optimization = TRUE),
        25.75,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'N1 works'
