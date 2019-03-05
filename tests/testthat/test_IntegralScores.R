context("Test integral scores")

test_that("Expected sample size is computed correctly",{
    null <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    dist <- Normal()

    # Define design from rpact
    design_rp2 <- rpact::getDesignInverseNormal(
        kMax = 2, alpha = 0.05, beta = 0.2, futilityBounds = 0, typeOfDesign = "P"
    )

    res    <- rpact::getSampleSizeMeans(
        design_rp2, normalApproximation = TRUE, alternative = .4
    )

    char   <- rpact::getDesignCharacteristics(design_rp2)

    n1 <- res$numberOfPatientsGroup1[1,]

    n2 <- res$numberOfPatientsGroup1[2,]

    c1f <- qnorm(char$futilityProbabilities) +
        sqrt(res$numberOfPatientsGroup1[1]) * .4 / sqrt(2)

    c1e <- design_rp2$criticalValues[1]

    f <- function(z){
        w1 <- 1 / sqrt(2)
        w2 <- sqrt(1 - w1^2)
        out <- (design_rp2$criticalValues[2] - w1 * z) / w2
        return(out)
    }

    x <- GaussLegendreRule(5)$nodes
    h <- (c1e - c1f) / 2
    x <- h * x + (h + c1f)

    design_gs <<- TwoStageDesign(
        round(n1),
        c1f,
        c1e,
        rep(round(n2), 5),
        sapply(seq(c1f, c1e, length.out = 5), f)
    )


    # Simulation
    sim_alt  <<- simulate(
        design_gs, nsim = 10000, dist = Normal(), theta = .4, seed = 59
    )
    sim_null <<- simulate(
        design_gs, nsim = 10000, dist = Normal(), theta = .0, seed = 59
    )


    # Expected Sample sizes under H1
    expect_equal(
        res$expectedPatientsH1 / 2, # per group!
        evaluate(expected(ConditionalSampleSize(dist, alternative)), design_gs),
        tolerance = 0.5
    ) # Compare ESS between specific package evaluation and rpact-value

    expect_equal(
        res$expectedPatientsH1 / 2,
        evaluate(expected(ConditionalSampleSize(dist, alternative)),
                 design_gs, specific = FALSE),
        tolerance = 0.5
    ) # compare ess between non-specific package evaluation and rpact-value

    expect_equal(
        res$expectedPatientsH1 / 2,
        mean(sim_alt[, "n1"]) + mean(sim_alt[,"n2"]),
        tolerance = 0.5
    ) # compare ESS between simluation result and rpact-value


    # Expected Sample sizes under H0
    expect_equal(
        res$expectedPatientsH0 / 2, # per group!
        evaluate(expected(ConditionalSampleSize(dist, null)), design_gs),
        tolerance = 0.5
    ) # compare ESS between specific package evaluation and rpact-value

    expect_equal(
        res$expectedPatientsH0 / 2,
        evaluate(expected(ConditionalSampleSize(dist, null)),
                 design_gs, specific = FALSE),
        tolerance = 0.5
    ) # compare ESS between non-specific package evaluation and rpact-value

    expect_equal(
        res$expectedPatientsH0 / 2,
        mean(sim_null[, "n1"]) + mean(sim_null[,"n2"]),
        tolerance = 0.5
    ) # compare ESS between simluation result and rpact-value

}) # end 'expected sample size is computed correctly'


test_that("Power is computed correctly for example design", {
    # Power
    pow <- expected(ConditionalPower(Normal(), PointMassPrior(.4, 1)))

    expect_equal(
        evaluate(pow, design_gs),
        .8,
        tolerance = .01
    )

    expect_equal(
        evaluate(pow, design_gs, specific = FALSE),
        .8,
        tolerance = .01
    )

    expect_equal(
        mean(sim_alt[, "reject"]),
        .8,
        tolerance = .01
    )
}) # end 'power is computed correctly for example design'


test_that("Type one error is computed correctly for example design", {
    # Type one error
    toer <- expected(ConditionalPower(Normal(), PointMassPrior(.0, 1)))

    expect_equal(
        evaluate(toer, design_gs),
        .05,
        tolerance = .005
    )

    expect_equal(
        evaluate(toer, design_gs, specific = FALSE),
        .05,
        tolerance = .005
    )

    expect_equal(
        mean(sim_null[, "reject"]),
        .05,
        tolerance = .005
    )

}) # end 'type one error is computed correctly for example design'



test_that("arithmetic works", {
    pow <- expected(ConditionalPower(Normal(), PointMassPrior(.35, 1)))
    ess <- expected(ConditionalSampleSize(Normal(), PointMassPrior(.35, 1)))

    expect_equal(
        evaluate(pow + ess, design_gs),
        evaluate(pow, design_gs) + evaluate(ess, design_gs)
    )

    expect_equal(
        evaluate(pow, design_gs) + 3,
        evaluate(pow + 3, design_gs)
    )

    expect_equal(
        1.1 + evaluate(pow, design_gs),
        evaluate(1.1 + pow, design_gs)
    )
}) # end 'arithmetic works'



test_that("show method returns class name", {
    pow <- expected(ConditionalPower(Normal(), PointMassPrior(.35, 1)))
    ess <- expected(ConditionalSampleSize(Normal(), PointMassPrior(.35, 1)))

    expect_equal(
        show(pow),
        show(ess)
    ) # should both be of class IntegralScore

}) # end 'show method returns class name'
