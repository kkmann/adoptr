context("TwoStageDesign")

test_that("Optimal design with point prior is computable", {
    # define an initial design
    n1     <- 25
    c1f    <-   .0
    c1e    <-  2.0
    number_knots <- 5
    n2_piv <- rep(40.0, number_knots)
    c2_piv <- rep(1.96, number_knots)
    design <- gq_design(n1, c1f, c1e, n2_piv, c2_piv, number_knots)

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
        2 * number_knots + 3
    )

    # check if key figures can be computed
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)

    dist <- Normal()

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

    #compute optimal design

    objective <- function(x) {
        d  <- update(design, x)
        evaluate(ess, d) + .001*evaluate(smth, d)
    }
    update(design, tunable_parameters(design))

    constraint <- function(x) {
        d  <- update(design, x)
        c(
            .8 - evaluate(pow, d),
            evaluate(toer, d) - 0.05,
            x[2] - x[3] + .1,
            diff(c2(d, scaled_integration_pivots(d)))
        )
    }

    ub <- c(50, 1, 4, numeric(number_knots) + 50, numeric(number_knots) + 5)
    lb <- c(10, -1, 1, numeric(number_knots) + 2, numeric(number_knots) - 5)

    res <- nloptr::nloptr(
        x0 = tunable_parameters(design),
        lb = lb,
        ub = ub,
        eval_f      = objective,
        eval_g_ineq = constraint,
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-4,
            maxeval     = 2500
        )
    )

    d2 <- update(design, res$solution)

    expect_equal(
        round(evaluate(pow, d2), 1),
        0.8
    )

    expect_equal(
        round(evaluate(toer, d2), 2),
        0.05
    )

    expect_equal(
        sign(evaluate(ess, d2) - evaluate(ess, design)),
        -1
    )


})
