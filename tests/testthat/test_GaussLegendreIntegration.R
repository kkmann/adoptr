context("Gauss Legendre Implementation")

test_that("Integration yields correct results", {

    # First integral, start with an easy one
    f1 <- function(x) .4 * x^5 + x^3 + 2*x + 1
    int1 <- stats::integrate(f1, -10, 20)$value
    rule1 <- GaussLegendreRule(3)
    gl1 <- GaussLegendreIntegral(f1, -10, 20, rule1)

    expect_equal(int1, gl1)


    # Let's become more difficult
    f2 <- function(x) sin(x) * cos(x)
    int2 <- stats::integrate(f2, -pi, pi/2)$value
    rule2 <- GaussLegendreRule(10)
    gl2 <- GaussLegendreIntegral(f2, -pi, pi/2, rule2)

    expect_equal(int2, gl2)


    # And now something that fits to our problem
    f3 <- function(x) pt(2^(3/2) - x, df = 30 - 2*x) * dt(x, df = 30 - 2*x)
    int3 <- stats::integrate(f3, -1, 3)$value
    rule3 <- GaussLegendreRule(10)
    gl3 <- GaussLegendreIntegral(f3, -1, 3, rule3)

    expect_equal(int3, gl3)


})
