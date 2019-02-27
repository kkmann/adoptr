context("tunable parameters")

# preliminaries
n1     <- 25
c1f    <-   .0
c1e    <-  2.0
order  <- 5L
n2_piv <- 40.0
c2_piv <- 1.96

initial_design <- TwoStageDesign(n1, c1f, c1e, n2_piv, c2_piv, order)

dist        <- Normal(two_armed = FALSE)
null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)

ess  <- expected(ConditionalSampleSize(dist, alternative))
pow  <- expected(ConditionalPower(dist, alternative))
toer <- expected(ConditionalPower(dist, null))





test_that("tunability of parameters can be changed", {

    tmp <- make_fixed(initial_design, n1, c1f, c1e)

    expect_true(!tmp@tunable["n1"])
    expect_true(!tmp@tunable["c1f"])
    expect_true(!tmp@tunable["c1e"])
    expect_equal(tmp@tunable["n2_pivots"], initial_design@tunable["n2_pivots"])
    expect_equal(tmp@tunable["c2_pivots"], initial_design@tunable["c2_pivots"])
    expect_equal(tmp@tunable["x1_norm_pivots"], initial_design@tunable["x1_norm_pivots"])
    expect_equal(tmp@tunable["weights"], initial_design@tunable["weights"])
    expect_equal(tmp@tunable["tunable"], initial_design@tunable["tunable"])

})



test_that("make_tunable works as desired",{
    tmp <- make_fixed(initial_design, n1, c1f, c1e, n2_pivots, c2_pivots)

    expect_true(!tmp@tunable["n1"])
    expect_true(!tmp@tunable["c1f"])
    expect_true(!tmp@tunable["c1e"])
    expect_true(!tmp@tunable["n2_pivots"])
    expect_true(!tmp@tunable["c2_pivots"])

    tmp <- make_tunable(initial_design, n1, c2_pivots)

    expect_true(tmp@tunable["n1"])
    expect_true(tmp@tunable["c2_pivots"])
    expect_equal(tmp@tunable["c2_pivots"], initial_design@tunable["c2_pivots"])

}) # end 'make_tunable works as desired'



test_that("two stage design can be optimized with fixed first stage", {

    tmp <- make_fixed(initial_design, n1, c1f, c1e)

    lb_design      <- update(tmp, c(numeric(order) + 1, numeric(order) - 3))
    ub_design      <- update(tmp, c(numeric(order) + 100, numeric(order) + 10))

    res <- suppressWarnings(minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design        = tmp,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design,
        opts = list(
            algorithm = "NLOPT_LN_COBYLA",
            maxiter   = 10 # we do not need convergence,
                           # only see if it works technically!
        )
    ))

    # check that fixed params did not change
    expect_equal(res$design@n1, tmp@n1)
    expect_equal(res$design@c1f, tmp@c1f)
    expect_equal(res$design@c1e, tmp@c1e)

})



test_that("two stage design can be optimized with fixed sample sizes", {

    tmp <- make_fixed(initial_design, n1, n2_pivots)

    lb_design      <- update(tmp, c(-1, 1, numeric(order) - 3))
    ub_design      <- update(tmp, c(1, 4, numeric(order) + 5))

    res <- suppressWarnings(minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design        = tmp,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design,
        opts = list(
            algorithm = "NLOPT_LN_COBYLA",
            maxiter   = 10 # we do not need convergence,
                           # only see if it works technically!
        )
    ))

    # check that fixed params did not change
    expect_equal(res$design@n1, tmp@n1)
    expect_equal(res$design@n2_pivots, tmp@n2_pivots)

})  # end 'two-stage design can be optimized with fixed sample sizes'




test_that("group-sequential design can be optimized with fixed sample sizes", {

    initial_gs_design <- GroupSequentialDesign(25, .0, 2.0, 50.0, 2.0, order)

    tmp_gs <- make_fixed(initial_gs_design, n1, n2_pivots)

    lb_design      <- update(tmp_gs, c(-1, 1, numeric(order) - 3))
    ub_design      <- update(tmp_gs, c(1, 4, numeric(order) + 5))

    res <- minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design        = tmp_gs,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design,
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_abs    = 1 # we do not need convergence,
                            # only see if it works technically!
        )
    )

    # make sure that the lax convergence still leads to function evaluations!
    expect_true(res$nloptr_return$iterations >= 10)

    # check that fixed params did not change
    expect_equal(res$design@n1, tmp_gs@n1)
    expect_equal(res$design@n2_pivots, tmp_gs@n2_pivots)

}) # end 'group-sequential design can be optimized with fixed sample sizes'



test_that("one-stage design can be optimized with fixed sample sizes", {

    initial_os_design <- OneStageDesign(60.0, 2.0)

    tmp_os <- make_fixed(initial_os_design, c1f)

    lb_design      <- update(tmp_os, 2.0)
    ub_design      <- update(tmp_os, 100.0)

    res <- minimize(

        ess,
        subject_to(
            pow >= 0.8,
            toer <= 0.025
            ),
        initial_design        = tmp_os,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design
    )

    # check that fixed params did not change
    expect_equal(res$design@c1f, tmp_os@c1f)

}) # end 'one-stage design can be optimized with fixed sample sizes'
