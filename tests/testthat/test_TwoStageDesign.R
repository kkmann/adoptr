context("TwoStageDesign")

test_that("gaussian quadrature constructor", {

    n1           <<-  49.6
    c1f          <<-   0.7
    c1e          <<-   2.5
    number_knots <<-   5L
    n2_piv       <<- rep(49.6, number_knots)
    c2_piv       <<- rep(1.96, number_knots)
    design       <<- TwoStageDesign(n1, c1f, c1e, n2_piv, c2_piv, number_knots)

    expect_equal(
        c(n1, c1f, c1e),
        c(design@n1, design@c1f, design@c1e),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    x1 <- seq(c1f, c1e, length.out = 11)

    expect_equal(
        n2(design, x1, round = FALSE),
        rep(49.6, length(x1)),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        n2(design, x1),
        rep(50, length(x1)),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        n(design, x1),
        rep(100, length(x1)),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'gaussian quadrature constructor'



test_that("simulate works (as last time)", {

    design@n1      <- 50

    expect_known_value(
        adoptr::simulate(design, nsim = 50, dist = Normal(), theta = .5, seed = 42),
        file = "known_values/simulate.rds")

}) # end 'simulate works'



test_that("errors are returned correctly", {

    # pivots length must fit
    expect_error(
        TwoStageDesign(50, 0, 2, rep(50, 3), c(2, 2)))

    cp      <- ConditionalPower(Normal(), PointMassPrior(.4, 1))
    pow     <- Power(Normal(), PointMassPrior(.4, 1))
    order   <- 5L
    design  <- TwoStageDesign(50.1, 0, 2, rep(50, order), rep(2, order))

    # unconditional scores are not plotted
    expect_error(
        plot(design, rounded = TRUE, "Power" = pow))

    # Conditional scores cannot be summarized
    expect_error(
        summary(design, rounded = TRUE, "Conditional Power" = cp))

}) # end 'errors are returned correctly'



test_that("print methods", {

    pow <- Power(Normal(), PointMassPrior(.4, 1))
    vdiffr::expect_doppelganger(
        "Design print",
        print.TwoStageDesignSummary(summary(design, "Power" = pow)))

}) # end 'print methods'



test_that("plot produces correct number of columns", {

    cp  <- ConditionalPower(Normal(), PointMassPrior(.3, 1))
    pic <- plot(design, "ConditionalPower" = cp)

    expect_true(
        pic$mfrow[2] == 3)

}) # end 'plot produces correct number of columns'



test_that("show method returns design name", {

    expect_true(
        class(design)[1] == capture.output(show(design)))

}) # end 'show method returns design name'



test_that("defining order does not destroy pivots", {

    n2 <- seq(100, 40, length.out = number_knots)
    c2 <- seq(2.0, 0.0, length.out = number_knots)
    d  <- TwoStageDesign(n1, c1f, c1e, n2, c2, number_knots)

    expect_equal(
        d@n2_pivots,
        n2,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        d@c2_pivots,
        c2,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'defining order does not destroy pivots'



test_that("boundary designs keep monotonicity", {

    n2   <- seq(100, 40, length.out = number_knots)
    c2   <- seq(2.0, 0.0, length.out = number_knots)
    d    <- TwoStageDesign(n1, c1f, c1e, n2, c2, number_knots)
    d_lb <- get_lower_boundary_design(d)
    d_ub <- get_upper_boundary_design(d)

    expect_true(all(
        sign(diff(d_lb@c2_pivots)) == sign(diff(d@c2_pivots))))

    expect_true(all(
        sign(diff(d_ub@n2_pivots)) == sign(diff(d@n2_pivots))))

    expect_true(all(
        sign(diff(d_ub@c2_pivots)) == sign(diff(d@c2_pivots))))

}) # end 'boundary designs keep monotonicity'
