context("FivePointDesign - optimization okay?                                ")

test_that("single point prior", {

    design <- FivePointDesign(25, 0, 2, rep(40.0, 5), rep( 1.96, 5))

    # define null and alternative as point mass distributions
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)

    dist <- Normal()

    ess  <- integrate(ConditionalSampleSize(dist, alternative))
    cp   <- ConditionalPower(dist, alternative)
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(dist, null))

    smth <- Smoothness_n2()

    objective <- function(x) {
        d  <- update(design, x)
        evaluate(ess, d) + .001*evaluate(smth, d)
    }

    constraint <- function(x) {
        d  <- update(design, x)
        c(
            .8 - evaluate(pow, d),
            evaluate(toer, d) - 0.05,
            x[2] - x[3] + .1,
            diff(c2(d, get_knots(d)))
        )
    }

    ub <- c(50, 1, 4, numeric(5) + 50, numeric(5) + 5)
    lb <- c(10, -1, 1, numeric(5) + 2, numeric(5) - 5)

    res <- nloptr::nloptr(
        as.numeric(design),
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

    expect_true(all(constraint(res$solution) < .0001))

})
