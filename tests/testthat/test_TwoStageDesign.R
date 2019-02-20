context("TwoStageDesign")

test_that("gaussian quadrature constructor", {

    n1           <<-  49.6
    c1f          <<-   0.7
    c1e          <<-   2.5
    number_knots <<-   5L
    n2_piv       <<- rep(49.6, number_knots)
    c2_piv       <<- rep(1.96, number_knots)
    design       <<- gq_design(n1, c1f, c1e, n2_piv, c2_piv, number_knots)

    expect_equal(
        c(n1, c1f, c1e),
        c(design@n1, design@c1f, design@c1e)
    )

    x1 <- seq(c1f, c1e, length.out = 11)

    expect_equal(
        n2(design, x1),
        rep(49.6, length(x1))
    )

    design@rounded <- TRUE # TODO: don't like this mechanism, why not make it an argument to n2()?

    expect_equal(
        n2(design, x1),
        rep(50, length(x1))
    )

    expect_equal(
        n(design, x1),
        rep(100, length(x1))
    )

}) # end 'gaussian quadrature constructor'



test_that("simulate works (as last time)", {

    design@n1      <- 50
    design@rounded <- TRUE

    expect_known_value(
        adoptr::simulate(design, nsim = 50, dist = Normal(), theta = .5, seed = 42),
        file = "known_values/simulate.rds"
    )

}) # end 'simulate works'



test_that("errors are returned correctly", {
    expect_error(
        TwoStageDesign(50, 0, 2, rep(50, 3), c(2, 2), c(.1, .5, .9), rep(1/3, 3))
    ) # pivots length must fit

    expect_error(
        TwoStageDesign(50, 0, 2, rep(50, 3), rep(2, 3), c(0, 1, 2), rep(1/3, 3))
    ) # x1_norm_pivots must not be scaled

    expect_error(
        TwoStageDesign(50, 0, 2, rep(50, 3), rep(2,3), c(-.5, 0, .5), c(1, 0, 1))
    ) # positive weights

    cp  <- ConditionalPower(Normal(), PointMassPrior(.4, 1))
    pow <- integrate(cp)
    order  = 5L
    design  <- gq_design(50.1, 0, 2, rep(50, order), rep(2, order), order)
    design2 <- gq_design(50, 0, 2, rep(50.1, order), rep(2, order), order)

    expect_error(
        plot(design, rounded = TRUE, "Power" = pow)
    ) # unconditional scores are not plotted

    expect_error(
        summary(design, rounded = TRUE, "Conditional Power" = cp)
    ) # Conditional scores cannot be summarized


    expect_error(
        simulate(design, nsim = 50, dist = Normal(), theta = .5, seed = 42)
    ) # simulate requires integer n1 values

    expect_error(
        simulate(design2, nsim = 50, dist = Normal(), theta = .5, seed = 42)
    ) # simulate requires integer n2 values

}) # end 'errors are returned correctly'
