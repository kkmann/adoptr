context("test point alternative point ESS design (case 01)")

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

    # what about increasing runtime?
    optimal_design <<- minimize(
        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design        = design,
        lower_boundary_design = update(design, c(10, -1, 1, numeric(order) + 2, numeric(order) - 5)),
        upper_boundary_design = update(design, c(50, 1, 4, numeric(order) + 50, numeric(order) + 5)),
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-5,
            maxeval     = 25000
        )
    )

    expect_known_value(optimal_design, file = "known_values/optimal_design_case_01.rds")

}) # end 'plausible n2() function'



