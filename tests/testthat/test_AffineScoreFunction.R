context("AffineScoreFunction")

# Define conditional and unconditional scores
css <- ConditionalSampleSize(Normal(), PointMassPrior(.3, 1))
cp  <- ConditionalPower(Normal(), PointMassPrior(.3, 1))
ess <- integrate(css)
pow <- integrate(cp)
order  = 5L
design <- gq_design(30, -.5, 2.5, rep(50.0, order), rep(2.0, order), order)
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
        afunsc + incpt,
        incpt + afunsc
    ) # addition is commutative

    expect_equal(
        afunsc * incpt,
        afunsc
    ) # multiplication by 1 does not change

    expect_equal(
        incpt * afunsc,
        afunsc
    ) # multipilcation is commutative

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
        afcosc + incpt,
        incpt + afcosc
    ) # addition is commutative

    expect_equal(
        afcosc * incpt,
        afcosc
    ) # multiplication by 1 does not change

    expect_equal(
        incpt * afcosc,
        afcosc
    ) # multipilcation is commutative

})

