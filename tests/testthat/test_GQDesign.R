context("GQDesign")

test_that("Optimal design with point prior is computable", {
    # define an initial design
    n1     <- 25
    c1f    <-   .0
    c1e    <-  2.0
    number_knots <- 5
    n2_piv <- rep(40.0, number_knots)
    c2_piv <- rep(1.96, number_knots)
    design <- GQDesign(n1, c1f, c1e, n2_piv, c2_piv, number_knots)

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
        length(as.numeric(design)),
        2 * number_knots + 3
    )

    # check if key figures can be computed
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)

    ess  <- integrate(SampleSize(alternative))
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(null))


    expect_equal(
        round(eval(ess, design), 1),
        44.1
    )

    expect_equal(
        round(eval(pow, design), 3),
        0.842
    )

    expect_equal(
        round(eval(toer, design), 3),
        0.035
    )

    #compute optimal design

    objective <- function(x) {
        d  <- update(design, x)
        eval(ess, d) + .001*eval(smth, d)
    }

    constraint <- function(x) {
        d  <- update(design, x)
        c(
            .8 - eval(pow, d),
            eval(toer, d) - 0.05,
            x[2] - x[3] + .1,
            diff(c2(d, get_knots(d)))
        )
    }

    ub <- c(50, 1, 4, numeric(number_knots) + 50, numeric(number_knots) + 5)
    lb <- c(10, -1, 1, numeric(number_knots) + 2, numeric(number_knots) - 5)

    res <- nloptr::nloptr(
        x0 = as.numeric(design),
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
        round(eval(pow, d2), 1),
        0.8
    )

    expect_equal(
        round(eval(toer, d2), 2),
        0.05
    )

    expect_equal(
        sign(eval(ess, d2) - eval(ess, design)),
        -1
    )


})
