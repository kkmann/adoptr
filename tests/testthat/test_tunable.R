context("tunable parameters")



n1     <- 25
c1f    <-   .0
c1e    <-  2.0
order  <- 5L
n2_piv <- rep(40.0, order)
c2_piv <- rep(1.96, order)
design <- gq_design(n1, c1f, c1e, n2_piv, c2_piv, order)

# check if key figures can be computed
null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)

dist <- Normal(two_armed = FALSE)

ess  <- integrate(ConditionalSampleSize(dist, alternative))
cp   <- ConditionalPower(dist, alternative)
pow  <- integrate(cp)
toer <- integrate(ConditionalPower(dist, null))
smth <- AverageN2()



test_that("tunability of parameters can be changed", {

    tmp <- make_fixed(design, n1, c1f, c1e)

    expect_true(!tmp@tunable["n1"])
    expect_true(!tmp@tunable["c1f"])
    expect_true(!tmp@tunable["c1e"])
    expect_equal(tmp@tunable["n2_pivots"], design@tunable["n2_pivots"])
    expect_equal(tmp@tunable["c2_pivots"], design@tunable["c2_pivots"])
    expect_equal(tmp@tunable["x1_norm_pivots"], design@tunable["x1_norm_pivots"])
    expect_equal(tmp@tunable["weights"], design@tunable["weights"])
    expect_equal(tmp@tunable["tunable"], design@tunable["tunable"])

})



test_that("two stage design can be optimized with fixed first stage", {

    tmp <- make_fixed(design, n1, c1f, c1e)


    suppressWarnings(
    opt_design <- minimize(
        ess + 0.0001 * smth,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design        = tmp,
        lower_boundary_design = update(tmp, c(numeric(order) + 2, numeric(order) - 5)),
        upper_boundary_design = update(tmp, c(numeric(order) + 50, numeric(order) + 5))
    )$design
    )

    expect_equal(opt_design@n1, tmp@n1)
    expect_equal(opt_design@c1f, tmp@c1f)
    expect_equal(opt_design@c1e, tmp@c1e)

})



test_that("two stage design can be optimized with fixed sample sizes", {

    # minimize design over all parameters
    suppressWarnings(
    opt_design <- minimize(
        ess + 0.0001*smth,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design = design,
        lower_boundary_design = update(design, c(10, -1, 1, numeric(order) + 2, numeric(order) - 5)),
        upper_boundary_design = update(design, c(50, 1, 4, numeric(order) + 50, numeric(order) + 5))
    )$design
    )

    # round sample sizes to nearest integers
    # TODO: we need an evalution mode where n2 is also rounded to the next integer!
    opt_design@n1        <- round(opt_design@n1)
    opt_design@n2_pivots <- round(opt_design@n2_pivots)

    tmp <- make_fixed(opt_design, n1, n2_pivots)

    # suppress initial design infeasible
    suppressWarnings(
    tmp <- minimize(
        ess + 0.0001 * smth,
        subject_to(
            pow  >= 0.8,
            toer <= .05
        ),
        initial_design        = tmp,
        lower_boundary_design = update(tmp, c(-1, 1, numeric(order) - 5)),
        upper_boundary_design = update(tmp, c(1, 4, numeric(order) + 5))
    )$design
    )

    expect_equal(opt_design@n1, tmp@n1)
    expect_equal(opt_design@n2_pivots, tmp@n2_pivots)

})



# TODO: Test that this also works for GS and OS designs!
