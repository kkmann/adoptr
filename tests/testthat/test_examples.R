context("test point alternative")

test_that("plausible n2() function", {

    order <- 9L

    design <- gq_design(
        n1 = 25,
        c1f = 0,
        c1e = 2,
        n2 = rep(40.0, order),
        c2 = rep( 1.96, order),
        order = order
    )

    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    datadist    <- Normal(two_armed = FALSE)

    ess  <- integrate(ConditionalSampleSize(datadist, alternative))
    cp   <- ConditionalPower(datadist, alternative)
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(datadist, null))

    optimal_design <- minimize(
        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design = design,
        lower_boundary_design = update(design, c(10, -1, 1, numeric(order) + 2, numeric(order) - 5)),
        upper_boundary_design = update(design, c(50, 1, 4, numeric(order) + 50, numeric(order) + 5)),
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-5,
            maxeval     = 2500, # increased precision
            print_level = 1
        )
    )

    plot(optimal_design)

    z1 <- seq(0, 2.5, by = .01)
    plot(
        z1,
        n(optimal_design, z1), ylim = c(0, 1.1*(optimal_design@n1 + max(optimal_design@n2_pivots))),
        'l'
    )
    lines(z1, 75*probability_density_function(datadist, z1, optimal_design@n1, .4))

    # okay, this does not make any sense at all, the alternative distribution
    # peaks at exactly the point where n() gets weird - something is going terribly
    # wrong here ;)

    expect_equal(1, 0) # just indicate that we have a problem ;)

    evaluate(ess, optimal_design)

    # manally fix the problem
    optimal_design@n2_pivots[order] <- 18

    evaluate(ess, optimal_design) # better result -> expected sample size could be ok

    # what about increasing runtime?
    optimal_design <- minimize(
        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design = optimal_design,
        lower_boundary_design = update(design, c(10, -1, 1, numeric(order) + 2, numeric(order) - 5)),
        upper_boundary_design = update(design, c(50, 1, 4, numeric(order) + 50, numeric(order) + 5)),
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-5,
            maxeval     = 10000, # increased precision
            print_level = 1
        )
    )

    plot(optimal_design)

    z1 <- seq(0, 2.5, by = .01)
    plot(
        z1,
        n(optimal_design, z1), ylim = c(0, 1.1*(optimal_design@n1 + max(optimal_design@n2_pivots))),
        'l'
    )

    # much better!


})
