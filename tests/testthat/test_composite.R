context("Score composition")



css    <- ConditionalSampleSize()
cp     <- ConditionalPower(Normal(), PointMassPrior(.3, 1))
ce     <- ConditionalPower(Normal(), PointMassPrior(.0, 1))
ess    <- expected(css, Normal(), PointMassPrior(.3, 1))
pow    <- expected(cp, Normal(), PointMassPrior(.3, 1))
order  <- 5L
design <- TwoStageDesign(30, -.5, 2.5, rep(50.0, order), rep(2.0, order))
z1     <- seq(-1, 3, .1)



test_that("Errors are defined correctly", {

    expect_error(
        composite({ess + cp}))

    expect_error(
        composite({1}))

}) # end 'errors are defined correctly'



test_that("Composition of unconditional scores", {

    expect_equal(
        evaluate(composite({2*ess + ess + 1}), design),
        3*evaluate(ess, design) + 1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(composite({ess + 10*pow}), design),
        evaluate(ess, design) + 10*evaluate(pow, design),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(composite({sin(ess)}), design),
        sin(evaluate(ess, design)),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'composition of unconditional scores'



test_that("Composition of conditional scores", {

    expect_equal(
        evaluate(composite({2*css + css + 1}), design, z1),
        3*evaluate(css, design, z1) + 1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(composite({css + 10*cp}), design, z1),
        evaluate(css, design, z1) + 10*evaluate(cp, design, z1),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(composite({sin(css)}), design, z1),
        sin(evaluate(css, design, z1)),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'composition of conditional scores'



test_that("Integrals of compositions", {

    expect_equal(
        evaluate(
            expected(composite({2*css + css + 1}), Normal(), PointMassPrior(.3, 1)),
            design
        ),
        3*evaluate(ess, design) + 1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'integrals of compositions'



test_that("Nested compositions", {

    cs1 <- composite({2*css})
    cs2 <- composite({3*cs1})
    cs3 <- function(a) composite({a*cs1})

    expect_equal(
        6 * evaluate(css, design, .5),
        evaluate(cs2, design, .5),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        evaluate(cs3(3), design, .5),
        3*evaluate(cs1, design, .5),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'nested compositions'



test_that("show method returns class name", {

    expect_equal(
        capture.output(show(composite({2*ess}))),
        "CompositeUnconditionalScore")

    expect_equal(
        capture.output(show(composite({2*css}))),
        "CompositeConditionalScore")

}) # end 'show method returns class name'
