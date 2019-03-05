context("AffineScoreFunction")

# Define conditional and unconditional scores
css <- ConditionalSampleSize(Normal(), PointMassPrior(.3, 1))
cp  <- ConditionalPower(Normal(), PointMassPrior(.3, 1))
ess <- expected(css)
pow <- expected(cp)
order  = 5L
design <- TwoStageDesign(30, -.5, 2.5, rep(50.0, order), rep(2.0, order))
z1 <- seq(-1, 3, .1)



test_that("Errors are defined correctly", {
    expect_error(
        AffineScore(c(ess, pow), 1, 0)
    )

    expect_error(
        AffineScore(c(ess, pow), c(1, 1), c(0, 1))
    )

    expect_error(
        AffineScore(c(ess, pow), c(NA, 1), 0)
    )

    expect_error(
        AffineConditionalScore(c(ess, pow), c(1, 1), 1)
    )

    expect_error(
        AffineUnconditionalScore(c(cp, css), c(1, 1), 1)
    )
}) # end 'errors are defined correctly'



test_that("Arithmetic for unconditional scores works", {
    afpow  <- AffineUnconditionalScore(pow, 1, 0)

    expect_equal(
        evaluate(afpow, design),
        evaluate(pow, design)
    ) # shifting by 0 does not change anything

    incpt  <- 1
    afunsc <- AffineUnconditionalScore(c(afpow, ess), c(.7, 1.3), -.1)

    expect_equal(
        evaluate(afunsc, design),
        .7 * evaluate(afpow, design) + 1.3 * evaluate(ess, design) - .1
    ) # arithmetic works

    expect_equal(
        afunsc + incpt,
        incpt + afunsc
    ) # addition with numeric is commutative

    expect_equal(
        afunsc * incpt,
        afunsc
    ) # multiplication by 1 does not change anything

    expect_equal(
        incpt * afunsc,
        afunsc
    ) # multiplication is commutative

    expect_equal(
        evaluate(2 * afpow + afunsc, design),
        2 * evaluate(afpow, design) + evaluate(afunsc, design)
    ) # add two AffineUnconditionalScores

}) # end 'Arithmetic for unconditional scores works'


test_that("Arithmetic for conditional scores works", {
    afcp   <- AffineConditionalScore(cp, 1, 0)

    expect_equal(
        evaluate(afcp, design, z1),
        evaluate(cp, design, z1)
    ) # shifting by 0 does not change anything

    incpt  <- 1
    afcosc <- AffineConditionalScore(c(afcp, css), c(20, .1), 1.1)

    expect_equal(
        evaluate(afcosc, design, z1),
        20 * evaluate(afcp, design, z1) + .1 * evaluate(css, design, z1) + 1.1
    ) # arithmetic works

    expect_equal(
        evaluate(afcp + css, design, z1),
        evaluate(css + afcp, design, z1)
    ) # addition is commutative

    expect_equal(
        afcosc + incpt,
        incpt + afcosc
    ) # addition with numeric is commutative

    expect_equal(
        afcosc * incpt,
        afcosc
    ) # multiplication by 1 does not change anything

    expect_equal(
        incpt * afcosc,
        afcosc
    ) # multiplication is commutative

    expect_equal(
        evaluate(-1 * afcp + afcosc, design, z1),
        -1 * evaluate(afcp, design, z1) + evaluate(afcosc, design, z1)
    ) # add two AffineConditionalScores

})




test_that("show method returns class name", {
    afpow  <- AffineUnconditionalScore(pow, 1, 0)
    afcp   <- AffineConditionalScore(cp, 1, 0)

    expect_equal(
        show(afpow),
        show(afcp)
    ) # should both be of class AffineScore

}) # end 'show method returns class name'

