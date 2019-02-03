context("Continuous Prior")

test_that("Optimal Design with Gaussian prior is computable", {

    # Define data distribution and type one error
    dist <- Normal(two_armed = T)
    null = PointMassPrior(0, 1)
    toer = integrate(ConditionalPower(dist, null))

    # Define alternative density
    density_1   <- function(x) dnorm(x, mean = .4, sd = .2)
    prior       <- ContinuousPrior(density_1, c(-.7, 1.6))
    alternative <- condition(prior, c(0, 1.6))

    # Power is conditional power on alternative
    pow <- integrate(ConditionalPower(dist, alternative))

    # Optimize on the whole parameter space by using the entire prior without conditioning
    ess_total <- integrate(ConditionalSampleSize(dist, prior))

    # Define an initial design
    order = 9L
    design <- gq_design(
        n1  = 50,
        c1f = .0,
        c1e = 2.0,
        n2  = seq(200, 20, length.out = order),
        c2  = rep(1.96, order),
        order = order)

    # Define smoothness term
    smth <- integrate(SmoothnessN2(dist))

    # Optimize
    suppressWarnings( #suppress that initial design is infeasible
    minimize(
        objective = ess_total + 1e-6*smth,
        subject_to(
            pow  >= 0.8,
            toer <= 0.025
        ),
        initial_design = design,
        lower_boundary_design = update(design, c(10, -2, 2, numeric(order) + 2, numeric(order) - 1)),
        upper_boundary_design = update(design, c(500, 2, 5, numeric(order) + 500, numeric(order) + 3))
    ) ->
        optimal_design
    )

    # Check outcome
    out <- summary(optimal_design, "Power" = pow, "Type One Error" = toer)

    # Check constraints
    expect_equal(
        round(as.numeric(out$scores["Power"]), 2),
        0.8
    )

    expect_equal(
        round(as.numeric(out$scores["Type One Error"]), 3),
        0.025
    )

    # Known from example
    expect_equal(
        sign(diff(optimal_design@c2_pivots)),
        rep(-1, (order - 1))
    )

    # Known from example
    expect_equal(
        optimal_design@n1,
        71,
        tolerance = 2
    )

})
