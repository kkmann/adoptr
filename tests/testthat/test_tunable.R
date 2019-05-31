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

ess  <- ExpectedSampleSize(dist, alternative)
pow  <- Power(dist, alternative)
toer <- Power(dist, null)


test_that("boundary designs respect tunable", {

    n2   <- seq(100, 40, length.out = order)
    c2   <- seq(2.0, 0.0, length.out = order)
    d    <- TwoStageDesign(n1, c1f, c1e, n2, c2, order)
    d    <- make_fixed(d, n1)
    d_lb <- get_lower_boundary_design(d)
    d_ub <- get_upper_boundary_design(d)

    expect_true(all(
        d@tunable == d_lb@tunable))

    expect_true(all(
        d@tunable == d_ub@tunable))

}) # end 'boundary designs respect tunable'



test_that("fixing parameters works", {

    # can change tunability of relevant parameters
    tmp <<- make_fixed(initial_design, n1, c1f, c1e, n2_pivots, c2_pivots)
    for (name in c("n1", "c1f", "c1e", "n2_pivots", "c2_pivots")) {
        expect_true(
            !tmp@tunable[name])
    }

    # did not affect other fields
    for (name in c("x1_norm_pivots", "weights")) {
        expect_true(
            tmp@tunable[name] == initial_design@tunable[name])
    }

})



test_that("fixed params can be made tunable again",{

    # can we make params 'tunable' again?
    tmp <<- make_tunable(initial_design, n1, c2_pivots)
    for (name in c("n1", "c1f", "c1e")) {
        expect_true(
            tmp@tunable[name])
    }

}) # end 'make_tunable works as desired'



test_that("two stage design can be optimized with fixed first stage", {

    initial_design <- make_fixed(initial_design, n1, c1f, c1e)

    res <- suppressWarnings(

        minimize(
            ess,
            subject_to(
                pow  >= 0.8,
                toer <= .05
            ),
            initial_design,
            opts = list(
                xtol_abs  = 1e-5,
                algorithm = "NLOPT_LN_COBYLA",
                maxiter   = 10 # we do not need convergence,
                               # only see if it works technically!
            )
        )
    )

    # check that fixed params did not change
    expect_equal(
        res$design@n1,
        initial_design@n1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        res$design@c1f,
        initial_design@c1f,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        res$design@c1e,
        initial_design@c1e,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

})



test_that("two stage design can be optimized with fixed sample sizes", {

    initial_design <- make_fixed(initial_design, n1, n2_pivots)

    res <- suppressWarnings(

        minimize(
             ess,
            subject_to(
                pow  >= 0.8,
                toer <= .05
            ),
            initial_design,
            opts = list(
                xtol_abs  = 1e-5,
                algorithm = "NLOPT_LN_COBYLA",
                maxiter   = 10 # we do not need convergence,
                               # only see if it works technically!
            )
        )
    )

    # check that fixed params did not change
    expect_equal(
        res$design@n1,
        initial_design@n1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        res$design@n2_pivots,
        initial_design@n2_pivots,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

})  # end 'two-stage design can be optimized with fixed sample sizes'




test_that("group-sequential design can be optimized with fixed sample sizes", {

    initial_gs_design <- make_fixed(
        GroupSequentialDesign(25, .0, 2.0, 50.0, 2.0, order),
        n1, n2_pivots)

    res <- suppressWarnings(

        minimize(
            ess,
            subject_to(
                pow  >= 0.8,
                toer <= .05
            ),
            initial_gs_design,
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_abs    = 1e-5,
                maxiter     = 10
            )
        )
    )

    # check that fixed params did not change
    expect_equal(
        res$design@n1,
        initial_gs_design@n1,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        res$design@n2_pivots,
        initial_gs_design@n2_pivots,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'group-sequential design can be optimized with fixed sample sizes'



test_that("one-stage design can be optimized with fixed sample sizes", {

    initial_os_design <- make_fixed(OneStageDesign(60.0, 2.0), c1f)

    res <- minimize(
        ess,
        subject_to(
            pow >= 0.8,
            toer <= 0.025
            ),
        initial_os_design
    )

    # check that fixed params did not change
    expect_equal(
        res$design@c1f,
        initial_os_design@c1f,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'one-stage design can be optimized with fixed sample sizes'
