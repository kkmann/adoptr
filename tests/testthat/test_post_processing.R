context("post processing")

test_that("Post processing yields to integer sample sizes", {
    # define an initial design
    order <- 5L # must be integer!

    design <- gq_design(
        n1 = 100,
        c1f = 0,
        c1e = 2,
        n2 = 10 * log(seq(600.0, 100.0, length.out = order)), #use something not too smooth
        c2 = rep( 1.96, order),
        order = order)

    # Try a continuous prior
    dichte_null <- function(x) dnorm(x, mean = 0, sd = .05)
    dichte_alt  <- function(x) dnorm(x, mean = .3, sd = .1)
    null        <- ContinuousPrior(dichte_null, c(-5, 5))
    alternative <- ContinuousPrior(dichte_alt, c(-5, 5))

    dist <- Normal()

    # Define key figures
    ess  <- integrate(ConditionalSampleSize(dist, alternative))
    cp   <- ConditionalPower(dist, alternative)
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(dist, null))
    smth <- integrate(SmoothnessN2(dist))


    #compute optimal design
    suppressWarnings( # suppress that initial design is infeasible
    minimize(
        objective = ess + 0.00001*smth,
        subject_to(
            pow  >= 0.8,
            toer <= .025
        ),
        initial_design = design,
        lower_boundary_design = update(design, c(10, -2, 1, numeric(order) + 10, numeric(order) - 2)),
        upper_boundary_design = update(design, c(500, 2, 4, numeric(order) + 500, numeric(order) + 3)),
        c2_monotone = F,
        post_process = T,
        opts = list(algorithm   = "NLOPT_LN_COBYLA",
                    xtol_rel    = 1e-4,
                    maxeval     = 3000)
    ) ->
        optimal_design
    )

    # Test outcome by using the specific summary function
    out <- summary(optimal_design, "power" = pow, "toer" = toer)

    expect_equal(
        out$design@n1,
        round(out$design@n1)
    )

    expect_equal(
        out$design@n2_pivots,
        round(out$design@n2_pivots)
    )

    expect_equal(
        round(as.numeric(out$scores["power"]), 1),
        0.8
    )

    expect_equal(
        round(as.numeric(out$scores["toer"]), 3),
        0.025
    )

})
