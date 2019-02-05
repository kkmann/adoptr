context("Test integral scores")

test_that("Expected sample size is computed correctly",{
    null <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    dist <- Normal()

    # Define design from rpact
    design <- rpact::getDesignInverseNormal(kMax = 2, alpha = 0.05, beta = 0.2,
                                            futilityBounds = 0, typeOfDesign = "P")
    res <- rpact::getSampleSizeMeans(design, normalApproximation = T, alternative = .4)
    char <- rpact::getDesignCharacteristics(design)
    n1 <- res$numberOfPatientsGroup1[1,]
    n2 <- res$numberOfPatientsGroup1[2,]
    f <- function(z){
        w1 <- 1 / sqrt(2)
        w2 <- sqrt(1 - w1^2)
        out <- (design$criticalValues[2] - w1 * z) / w2
        return(out)
    }
    c1f <- qnorm(char$futilityProbabilities) + sqrt(res$numberOfPatientsGroup1[1]) * .4 / sqrt(2)
    c1e <- design$criticalValues[1]
    x <- GaussLegendreRule(5)$nodes
    h <- (c1e - c1f) / 2
    x <- h * x + (h + c1f)
    design_gs <- gq_design(round(n1),
                           c1f,
                           c1e,
                           rep(round(n2), 5),
                           sapply(seq(c1f, c1e, length.out = 5), f),
                           5L)


    # Simulation
    sim_alt <- simulate(design_gs, nsim = 10000, dist = Normal(), theta = .4, seed = 59)
    sim_null <- simulate(design_gs, nsim = 10000, dist = Normal(), theta = .0, seed = 59)


    # Expected Sample sizes under H1
    expect_equal(
        res$expectedPatientsH1 / 2, # per group!
        evaluate(integrate(ConditionalSampleSize(dist, alternative)), design_gs),
        tolerance = 0.5)

    expect_equal(res$expectedPatientsH1 / 2,
                 evaluate(integrate(ConditionalSampleSize(dist, alternative)), design_gs, specific = F),
                 tolerance = 0.5)

    expect_equal(res$expectedPatientsH1 / 2,
                 mean(sim_alt[, "n1"]) + mean(sim_alt[,"n2"]),
                 tolerance = 0.5)


    # Expected Sample sizes under H0
    expect_equal(
        res$expectedPatientsH0 / 2, # per group!
        evaluate(integrate(ConditionalSampleSize(dist, null)), design_gs),
        tolerance = 0.5)
    expect_equal(res$expectedPatientsH0 / 2,
                 evaluate(integrate(ConditionalSampleSize(dist, null)), design_gs, specific = F),
                 tolerance = 0.5)
    expect_equal(res$expectedPatientsH0 / 2,
                 mean(sim_null[, "n1"]) + mean(sim_null[,"n2"]),
                 tolerance = 0.5)

}) # end 'expected sample size is computed correctly'


test_that("Error rates are computed correctly", {
    # Power
    pow <- integrate(ConditionalPower(dist, alternative))
    expect_equal(
        evaluate(pow, design_gs),
        .8,
        tolerance = .01)
    expect_equal(
        evaluate(pow, design_gs, specific = F),
        .8,
        tolerance = .01)
    expect_equal(
        mean(sim_alt[, "reject"]),
        .8,
        tolerance = .01)


    # Type one error
    toer <- integrate(ConditionalPower(dist, null))
    expect_equal(
        evaluate(toer, design_gs),
        .05,
        tolerance = .005)
    expect_equal(
        evaluate(toer, design_gs, specific = F),
        .05,
        tolerance = .005)
    expect_equal(
        mean(sim_null[, "reject"]),
        .05,
        tolerance = .005)

}) # end 'error rates are computed correctly'
