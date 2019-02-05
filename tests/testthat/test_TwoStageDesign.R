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
        simulate(design, nsim = 50, dist = Normal(), theta = .5, seed = 42),
        file = "known_values/simulate.rds"
    )

}) # end 'simulate works'
