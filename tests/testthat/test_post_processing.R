context("post processing")

test_that("Post processing yields to integer sample sizes", {
    # define an initial design
    order <- 5L # must be integer!

    design <- gq_design(
        n1  = 100,
        c1f = 0,
        c1e = 2,
        n2  = 10 * log(seq(600.0, 100.0, length.out = order)), #use something not too smooth
        c2  = rep( 1.96, order),
        order = order)

    # Try a continuous prior
    null        <- PointMassPrior(.0, 1)
    dichte_alt  <- function(x) dnorm(x, mean = .4, sd = .1)
    alternative <- ContinuousPrior(dichte_alt, c(-.5, 1.5))
    alternative <- condition(alternative, c(0, 1.5))


    dist <- Normal(two_armed = FALSE)

    # Define key figures
    ess  <- integrate(ConditionalSampleSize(dist, alternative))
    cp   <- ConditionalPower(dist, alternative)
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(dist, null))
    avn2 <- AverageN2()


    #compute optimal design
    suppressWarnings( # suppress that initial design is infeasible
    minimize(
        objective = ess + 0.000001*avn2,
        subject_to(
            pow  >= 0.8,
            toer <= .025
        ),
        initial_design = design,
        lower_boundary_design = update(design, c(10, -1, 2, numeric(order) + 10, numeric(order) - 2)),
        upper_boundary_design = update(design, c(500, 2, 4, numeric(order) + 500, numeric(order) + 3)),
        c2_monotone  = FALSE,
        post_process = TRUE,
        opts = list(algorithm   = "NLOPT_LN_COBYLA",
                    xtol_rel    = 1e-4,
                    maxeval     = 3000)
    ) ->
        optimal_design
    )

    # Test outcome by using the specific summary function
    out <- summary(optimal_design$design, "power" = pow, "toer" = toer)

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

    expect_equal(
        optimal_design$design@c2_pivots,
        optimal_design$details$nloptr_return_post$solution[3 : (2+order)]
    ) # test if nloptr output works

})



test_that("nloptr warnings are returned correctly", {
    order <- 5L

    design <- gq_design(25, 0, 2, rep(40.0, order), rep(1.96, order), order)
    lb_design <- update(design, c(10, -1, 2, numeric(order) + 2, numeric(order) - 2))
    ub_design <- update(design, c(500, 2, 5, numeric(order) + 500, numeric(order) + 5))

    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    datadist    <- Normal(two_armed = TRUE)

    ess   <- integrate(ConditionalSampleSize(datadist, alternative))
    cp    <- ConditionalPower(datadist, alternative)
    pow   <- integrate(cp)
    toer  <- integrate(ConditionalPower(datadist, null))


    expect_warning(
        minimize(
            ess,
            subject_to(
                pow  >= 0.8,
                toer <= 0.05
            ),
            initial_design        = design,
            lower_boundary_design = lb_design,
            upper_boundary_design = ub_design,
            post_process = TRUE,
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_rel    = 1e-5,
                maxeval     = 100
            )
        )
    )

}) # end 'nloptr warnings are returned correctly'
