context("manual Gauss-Legendre integration")



test_that("manual Gauss-Legendre integration is replicating stats::integrate()", {

    f1    <- function(x) .4 * x^5 + x^3 + 2*x + 1
    int1  <- stats::integrate(f1, -10, 20)$value
    rule1 <- GaussLegendreRule(3)
    gl1   <- integrate_rule(f1, -10, 20, rule1$nodes, rule1$weights)

    expect_equal(
        int1,
        gl1,
        tolerance = 1e-6, scale = 1)

    f2    <- function(x) sin(x) * cos(x)
    int2  <- stats::integrate(f2, -pi, pi/2)$value
    rule2 <- GaussLegendreRule(10)
    gl2   <- integrate_rule(f2, -pi, pi/2, rule2$nodes, rule2$weights)

    expect_equal(
        int2,
        gl2,
        tolerance = 1e-6, scale = 1)

    f3 <- function(x) pt(2^(3/2) - x, df = 30 - 2*x) * dt(x, df = 30 - 2*x)
    int3 <- stats::integrate(f3, -1, 3)$value
    rule3 <- GaussLegendreRule(10)
    gl3 <- integrate_rule(f3, -1, 3, rule3$nodes, rule3$weights)

    expect_equal(
        int3,
        gl3,
        tolerance = 1e-6, scale = 1)

})



test_that("errors are thrown correctly", {

    expect_error(
        GaussLegendreRule(0))

    f <- function(x) x^2
    a <- 0
    b <- 1

    expect_error(
        integrate_rule(f, a, b, c(1, 1), 1))

    expect_error(
        integrate_rule(f, a, b, c(2, 1), c(.5, .5)))


    expect_error(
        integrate_rule(f, a, b, c(-.5, .5), c(-1, 1)))

}) # end 'errors are thrown correctly'
