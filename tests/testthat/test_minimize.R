context("check minimize()")



#preliminaries
order <- 5L

initial_design <- gq_design(25, 0, 2, rep(40.0, order), rep(1.96, order), order)
lb_design      <- update(design, c(5, -1, 2, numeric(order), numeric(order) - 3))
ub_design      <- update(design, c(100, 2, 5, numeric(order) + 100, numeric(order) + 5))

null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)
datadist    <- Normal(two_armed = TRUE)

ess   <- integrate(ConditionalSampleSize(datadist, alternative))
pow   <- integrate(ConditionalPower(datadist, alternative))
toer  <- integrate(ConditionalPower(datadist, null))



test_that("nloptr maxiter warning correctly passed", {

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
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_rel    = 1e-5,
                maxeval     = 10
            )
        )
    )
})



test_that("nloptr invalid initial values error works", {

    lb_design@n1 <- design@n1 + 1

    expect_error(
        minimize(
            ess,
            subject_to(
                pow  >= 0.8,
                toer <= 0.05
            ),
            initial_design        = design,
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
