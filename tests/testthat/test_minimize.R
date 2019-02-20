context("check minimize()")



#preliminaries
order <- 5L

initial_design <- gq_design(25, 0, 2, rep(40.0, order), rep(1.96, order), order)
lb_design      <- update(initial_design, c(5, -1, 2, numeric(order), numeric(order) - 3))
ub_design      <- update(initial_design, c(100, 2, 5, numeric(order) + 100, numeric(order) + 5))

null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)
datadist    <- Normal(two_armed = FALSE)

ess   <- integrate(ConditionalSampleSize(datadist, alternative))
cp    <- ConditionalPower(datadist, alternative)
pow   <- integrate(cp)
toer  <- integrate(ConditionalPower(datadist, null))



test_that("nloptr maxiter warning correctly passed", {

    expect_warning(
        minimize(
            ess,
            subject_to(
                pow  >= 0.8,
                toer <= 0.05
            ),

            post_process          = FALSE,
            initial_design        = initial_design,
            lower_boundary_design = lb_design,
            upper_boundary_design = ub_design,
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_rel    = 1e-5,
                maxeval     = 10
            )
        )
    )
})



test_that("nloptr invalid initial values error works", {

    lb_design@n1 <- initial_design@n1 + 1

    expect_error(
        minimize(
            ess,
            subject_to(
                pow  >= 0.8,
                toer <= 0.05
            ),

            post_process          = FALSE,
            initial_design        = initial_design,
            lower_boundary_design = lb_design,
            upper_boundary_design = ub_design,
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_rel    = 1e-5,
                maxeval     = 10
            )
        )
    )

})



test_that("post-processing yields integer sample sizes", {

    res <- suppressWarnings(minimize( # we do not warning about non-convergence!

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),

        initial_design        = initial_design,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design,
        post_process          = TRUE,
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-4, # we use drastically reduced precision, not about convergence!
            maxeval     = 10
        )

    ))

    # n1 is integer
    expect_equal(
        res$design@n1,
        round(res$design@n1)
    )

    # n2 is integer
    expect_equal(
        res$design@n2_pivots,
        round(res$design@n2_pivots)
    )

})



test_that("base-case satisfies constraints", {

    res <- minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),

        initial_design        = initial_design,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design,
        post_process          = TRUE,
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-3, # we use reduced precision, not optimal but should
                                # respect constraints!
            maxeval     = 1000
        )

    )

    # compute summaries
    out <- summary(res$design, "power" = pow, "toer" = toer)

    expect_equal(
        round(as.numeric(out$scores["power"]), 1),
        0.8
    )

    expect_equal(
        round(as.numeric(out$scores["toer"]), 3),
        0.05
    )

}) # end base-case respects constraints


test_that("base-case results are consistent", {

    opt_os <- minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),

        post_process          = FALSE,
        initial_design        = OneStageDesign(100, 1.97),
        lower_boundary_design = OneStageDesign(1, -5),
        upper_boundary_design = OneStageDesign(200, 5)

    )



    initial_design_gs <- gq_design(25, 0, 2, 40, rep(1.96, order), order)

    opt_gs <- minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),

        post_process          = TRUE,
        initial_design        = initial_design_gs,
        lower_boundary_design = update(initial_design_gs, c(10, -1, 1, 2, numeric(order) - 5)),
        upper_boundary_design = update(initial_design_gs, c(50, 1, 4, 50, numeric(order) + 5))

    )

    opt_ts <- minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),

        post_process          = TRUE,
        initial_design        = initial_design,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design

    )

    # TODO: check vs simulated results!

    # TODO: check for consistency ts better gs better os!

})


# TODO: one more case (no post processing) checking that conditional constraints work too!
