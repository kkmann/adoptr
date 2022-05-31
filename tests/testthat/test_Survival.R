context("Survival class")

test_that("Constructor works",  {

    expect_true(
        Survival(0.4)@two_armed)

    expect_true(
        !Survival(0.6,two_armed=FALSE)@two_armed)

    expect_error(
        Survival(-.1))

    expect_error(
        Survival(1.1))
})

test_that("pdf is defined correctly" , {

    dist <- Survival(0.7)
    x       <- seq(-3, 3, by = .1)
    n       <- 22
    theta   <- 1.3
    expect_equal(
        probability_density_function(dist,x,n,theta),
        stats::dnorm(x,log(theta)*sqrt(n/2),1),
        tolerance = 1e-6, scale = 1)
})

test_that("cdf is defined correctly" , {

    dist <- Survival(0.8)
    x       <- seq(-3, 3, by = .1)
    n       <- 11
    theta   <- 0.9
    expect_equal(
        cumulative_distribution_function(dist,x,n,theta),
        stats::pnorm(x,log(theta)*sqrt(n/2),1),
        tolerance = 1e-6, scale = 1)
})

test_that("quantile is defined correctly" , {

    dist <- Survival(0.4)
    x       <- seq(0.1, 1, by = .01)
    n       <- 44
    theta   <- 1.7
    expect_equal(
        quantile(dist,x,n,theta),
        stats::qnorm(x,log(theta)*sqrt(n/2),1),
        tolerance = 1e-6, scale = 1)

    expect_warning(
        quantile(dist, -1, n, theta)
    )
})

test_that("simulate respects seed", {

    expect_equal(
        simulate(Survival(0.7), 10, 10, 1, seed = 42),
        simulate(Survival(0.7), 10, 10, 1, seed = 42),
        tolerance = 1e-6, scale = 1)

    set.seed(42)

    expect_true(
        all(simulate(Survival(0.7), 10, 12, 1.1) != simulate(Survival(0.7), 10, 12, 1.1)))
})

test_that("Survival and Normal equal under null hypothesis", {
    dist1 <- Survival(0.7)
    dist2 <- Normal()
    x <- seq(-3,3,by=0.1)
    n <- 10
    expect_equal(
        probability_density_function(dist1,x,n,1),
        probability_density_function(dist2,x,n,0),
        tolerance = 1e-6, scale = 1
    )
})

test_that("show method", {

    expect_equal(
        capture.output(show(Survival(0.7))),
        "Survival<two-armed>, event rate: 0.7 "
    )

    expect_equal(
        capture.output(show(Survival(0.7, two_armed = FALSE))),
        "Survival<single-armed>, event rate: 0.7 "
    )

})
