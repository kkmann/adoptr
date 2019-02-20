context("tunable parameters")

# preliminaries
n1     <- 25
c1f    <-   .0
c1e    <-  2.0
order  <- 5L
n2_piv <- rep(40.0, order)
c2_piv <- rep(1.96, order)

initial_design <- gq_design(n1, c1f, c1e, n2_piv, c2_piv, order)

dist        <- Normal(two_armed = FALSE)
null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)

ess  <- integrate(ConditionalSampleSize(dist, alternative))
pow  <- integrate(ConditionalPower(dist, alternative))
toer <- integrate(ConditionalPower(dist, null))





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

# TODO: does make_tunable work as well?



test_that("two stage design can be optimized with fixed first stage", {

    tmp <- make_fixed(initial_design, n1, c1f, c1e)

    lb_design      <- update(tmp, c(numeric(order) + 1, numeric(order) - 3))
    ub_design      <- update(tmp, c(numeric(order) + 100, numeric(order) + 10))

    res <- minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),

        post_process          = FALSE, # if we can do it without post
                                       # processing, by design, it also
                                       # works with
        initial_design        = tmp,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design,
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_abs    = 1 # we do not need convergence,
                            # only see if it works technically!
        )
    )

    # make sure that the lax convergence still leads to function evaluations!
    expect_true(res$details$nloptr_return$iterations >= 10)

    # check that fixed params did not change
    expect_equal(res$design@n1, tmp@n1)
    expect_equal(res$design@c1f, tmp@c1f)
    expect_equal(res$design@c1e, tmp@c1e)

})



test_that("two stage design can be optimized with fixed sample sizes", {

    tmp <- make_fixed(initial_design, n1, n2_pivots)

    lb_design      <- update(tmp, c(-1, 1, numeric(order) - 3))
    ub_design      <- update(tmp, c(1, 4, numeric(order) + 5))

    res <- minimize(

        ess,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),

        post_process          = FALSE, # if we can do it without post
                                       # processing, by design, it also
                                       # works with
        initial_design        = tmp,
        lower_boundary_design = lb_design,
        upper_boundary_design = ub_design,
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_abs    = 1 # we do not need convergence,
                            # only see if it works technically!
        )

    )

    # make sure that the lax convergence still leads to function evaluations!
    expect_true(res$details$nloptr_return$iterations >= 10)

    # check that fixed params did not change
    expect_equal(res$design@n1, tmp@n1)
    expect_equal(res$design@n2_pivots, tmp@n2_pivots)

})



# TODO: Test that this also works for GS and OS designs!
