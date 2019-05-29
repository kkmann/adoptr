context("Test integral scores")



expect_delta_within <- function(x1, x2, abs_tol) expect_lt(abs(x1 - x2), abs_tol)

tol_n   <- .75 # tolerance for sample size differneces
tol_pow <- .01 # tolerance for power differences
tol_a   <- .001 # tolerance for type one error rates



test_that("Expected sample size is computed correctly",{

    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.4, 1)
    dist        <- Normal()

    # Define design from rpact
    design_rp2 <- rpact::getDesignInverseNormal(
        kMax = 2, alpha = 0.05, beta = 0.2, futilityBounds = 0, typeOfDesign = "P"
    )

    res  <- rpact::getSampleSizeMeans(
        design_rp2, normalApproximation = TRUE, alternative = .4
    )
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
        round(n1),
        c1f,
        c1e,
        rep(round(n2), 5),
        sapply(x, f)
    )


    # Simulation
    sim_alt  <<- simulate(
        design_gs, nsim = 10000, dist = Normal(), theta = .4, seed = 59
    )
    sim_null <<- simulate(
        design_gs, nsim = 10000, dist = Normal(), theta = .0, seed = 59
    )


    # Expected Sample sizes under H1
    expect_delta_within(
        res$expectedNumberOfSubjectsH1/2, # per group!
        evaluate(expected(ConditionalSampleSize(), dist, alternative), design_gs),
        tol_n) # Compare ESS between specific package evaluation and rpact-value

    expect_delta_within(
        res$expectedNumberOfSubjectsH1/2,
        evaluate(expected(ConditionalSampleSize(), dist, alternative),
                   design_gs, specific = FALSE),
        tol_n) # compare ess between non-specific package evaluation and rpact-value

    expect_delta_within(
        res$expectedNumberOfSubjectsH1/2,
        (mean(sim_alt[, "n1"]) + mean(sim_alt[,"n2"])),
        tol_n) # compare ESS between simluation result and rpact-value


    # Expected Sample sizes under H0
    expect_delta_within(
        res$expectedNumberOfSubjectsH1/2,
        evaluate(expected(ConditionalSampleSize(), dist, null), design_gs),
        tol_n) # compare ESS between specific package evaluation and rpact-value

    expect_delta_within(
        res$expectedNumberOfSubjectsH0/2,
        evaluate(expected(ConditionalSampleSize(), dist, null),
                 design_gs, specific = FALSE),
        tol_n) # compare ESS between non-specific package evaluation and rpact-value

    expect_delta_within(
        res$expectedNumberOfSubjectsH0/2,
        mean(sim_null[, "n1"]) + mean(sim_null[,"n2"]),
        tol_n) # compare ESS between simluation result and rpact-value

}) # end 'expected sample size is computed correctly'



test_that("Power is computed correctly for example design", {

    # Power
    pow <- Power(Normal(), PointMassPrior(.4, 1))

    expect_delta_within(evaluate(pow, design_gs), .8, tol_pow)
    expect_delta_within(evaluate(pow, design_gs, specific = FALSE), .8, tol_pow)
    expect_delta_within(mean(sim_alt[, "reject"]), .8, tol_pow)

}) # end 'power is computed correctly for example design'



test_that("Type one error is computed correctly for example design", {
    # Type one error
    toer <- Power(Normal(), PointMassPrior(.0, 1))

    expect_delta_within(evaluate(toer, design_gs), .05, tol_a)
    expect_delta_within(evaluate(toer, design_gs, specific = FALSE), .05, tol_a)
    expect_delta_within(mean(sim_null[, "reject"]), .05, tol_a)

}) # end 'type one error is computed correctly for example design'
