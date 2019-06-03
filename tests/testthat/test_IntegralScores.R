context("Test integral scores")



test_that("Expected sample size is computed correctly",{

    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    dist        <- Normal()

    # Define design from rpact
    design_rp2 <- rpact::getDesignInverseNormal(
        kMax = 2, alpha = 0.05, beta = 0.2, futilityBounds = 0, typeOfDesign = "P")

    res  <- rpact::getSampleSizeMeans(
        design_rp2, normalApproximation = TRUE, alternative = .4)
    char <- rpact::getDesignCharacteristics(design_rp2)
    n1   <- res$numberOfSubjects1[1, 1]
    n2   <- res$numberOfSubjects1[2, 1] - n1
    c1f  <- qnorm(char$futilityProbabilities) +
        sqrt(res$numberOfSubjects1[1]) * .4 / sqrt(2)
    c1e  <- design_rp2$criticalValues[1]

    f <- function(z){
        w1 <- 1 / sqrt(2)
        w2 <- sqrt(1 - w1^2)
        out <- (design_rp2$criticalValues[2] - w1 * z) / w2
        return(out)
    }

    x <- adoptr:::GaussLegendreRule(5)$nodes
    h <- (c1e - c1f) / 2
    x <- h * x + (h + c1f)

    design_gs <<- TwoStageDesign(
        n1,
        c1f,
        c1e,
        rep(n2, 5),
        sapply(x, f))


    # Simulation
    sim_alt  <<- simulate(
        design_gs, nsim = 1e4, dist = Normal(), theta = .4, seed = 59)

    sim_null <<- simulate(
        design_gs, nsim = 1e4, dist = Normal(), theta = .0, seed = 59)


    # optimization = TRUE uses non-rounded values of n (as does rpact!)

    # Expected Sample sizes under H1
    expect_equal(
        res$expectedNumberOfSubjectsH1/2, # per group!
        evaluate(ExpectedSampleSize(dist, alternative), design_gs, optimization = TRUE),
        tolerance = .1, scale = 1)

    expect_equal(
        res$expectedNumberOfSubjectsH1/2,
        evaluate(ExpectedSampleSize(dist, alternative),
                   design_gs, specific = FALSE, optimization = TRUE),
        tolerance = .1, scale = 1)

    # Expected Sample sizes under H0
    expect_equal(
        res$expectedNumberOfSubjectsH1/2,
        evaluate(expected(ConditionalSampleSize(), dist, null), design_gs, optimization = TRUE),
        tolerance = .1, scale = 1)

    expect_equal(
        res$expectedNumberOfSubjectsH0/2,
        evaluate(expected(ConditionalSampleSize(), dist, null),
                 design_gs, specific = FALSE, optimization = TRUE),
        tolerance = .1, scale = 1)

}) # end 'expected sample size is computed correctly'



test_that("Power is computed correctly for example design", {

    pow <- Power(Normal(), PointMassPrior(.4, 1))

    expect_equal(
        evaluate(pow, design_gs),
        .8, ,
        tolerance = 1e-2, scale = 1)

    expect_equal(
        evaluate(pow, design_gs, specific = FALSE),
        .8,
        tolerance = 1e-2, scale = 1)

    expect_equal(
        mean(sim_alt[, "reject"]),
        .8, ,
        tolerance = 1e-2, scale = 1)

}) # end 'power is computed correctly for example design'



test_that("Type one error is computed correctly for example design", {

    toer <- Power(Normal(), PointMassPrior(.0, 1))
    expect_equal(
        evaluate(toer, design_gs),
        .05,
        tolerance = 1e-3, scale = 1)

    expect_equal(
        evaluate(toer, design_gs, specific = FALSE),
        .05,
        tolerance = 1e-3, scale = 1)

    expect_equal(
        mean(sim_null[, "reject"]),
        .05,
        tolerance = 1e-3, scale = 1)

}) # end 'type one error is computed correctly for example design'
