context("Score composition")

# Define conditional and unconditional scores
css    <- ConditionalSampleSize(Normal(), PointMassPrior(.3, 1))
cp     <- ConditionalPower(Normal(), PointMassPrior(.3, 1))
ce     <- ConditionalPower(Normal(), PointMassPrior(.0, 1))
ess    <- expected(css)
pow    <- expected(cp)

order  <- 5L
design <- TwoStageDesign(30, -.5, 2.5, rep(50.0, order), rep(2.0, order))
z1     <- seq(-1, 3, .1)



test_that("Errors are defined correctly", {
    expect_error(
        compose({ess + cp})
    )

    expect_error(
        compose({1})
    )

}) # end 'errors are defined correctly'



test_that("Composition of unconditional scores", {

    expect_equal(
        evaluate(compose({2*ess + ess + 1}), design),
        3*evaluate(ess, design) + 1
    )

    expect_equal(
        evaluate(compose({ess + 10*pow}), design),
        evaluate(ess, design) + 10*evaluate(pow, design)
    )

    expect_equal(
        evaluate(compose({sin(ess)}), design),
        sin(evaluate(ess, design))
    ) # functional composition

}) # end 'Linear combinations of unconditional scores work'



test_that("Composition of conditional scores", {

    expect_equal(
        evaluate(compose({2*css + css + 1}), design, z1),
        3*evaluate(css, design, z1) + 1
    )

    expect_equal(
        evaluate(compose({css + 10*cp}), design, z1),
        evaluate(css, design, z1) + 10*evaluate(cp, design, z1)
    )

    expect_equal(
        evaluate(compose({sin(css)}), design, z1),
        sin(evaluate(css, design, z1))
    ) # functional composition

})



test_that("Integrals of compositions", {

    expect_equal(
        evaluate(expected(compose({2*css + css + 1})), design),
        3*evaluate(ess, design) + 1
    )

}) # end 'show method returns class name'




test_that("Nested compositions", {

    s <- compose({2*css + css + 1})

    expect_equal(
        evaluate(expected(compose({2*css + css + 1})), design),
        3*evaluate(ess, design) + 1
    )

}) # end 'show method returns class name'



test_that("show method returns class name", {

    expect_equal(
        capture.output(show(compose({2*ess}))),
        "CompositeUnconditionalScore"
    )

    expect_equal(
        capture.output(show(compose({2*css}))),
        "CompositeConditionalScore"
    )

}) # end 'show method returns class name'

