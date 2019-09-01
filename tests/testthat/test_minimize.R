context("check minimize()")



# preliminaries
order <- 5L

initial_design <- TwoStageDesign(25, 0, 2, rep(35.0, order), rep(1.96, order))

null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)
datadist    <- Normal(two_armed = FALSE)

ess   <- ExpectedSampleSize(datadist, alternative)
ess_0 <- ExpectedSampleSize(datadist, null)
cp    <- ConditionalPower(datadist, alternative)
pow   <- expected(cp, datadist, alternative)
toer  <- Power(datadist, null)

alpha <- 0.05
beta  <- 0.2



test_that("nloptr maxiter warning correctly passed", {

    expect_warning(

        minimize(
            ess,
            subject_to(
                pow  >= 1 - beta,
                toer <= alpha
            ),
            initial_design,
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_rel    = 1e-5,
                maxeval     = 10
            )
        )
    )

})



test_that("nloptr invalid initial values error works", {

    lb_design    <- update(initial_design, c(5, -1, 2, numeric(order), numeric(order) - 3))
    lb_design@n1 <- initial_design@n1 + 1

    expect_error(

        minimize(
            ess,
            subject_to(
                pow  >= 1 - beta,
                toer <= alpha
            ),
            initial_design,
            lower_boundary_design = lb_design,
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_rel    = 1e-5,
                maxeval     = 10
            )
        )

    )

})


test_that("Optimal one-stage design can be computed", {

    opt_os <<- minimize(
        ess,
        subject_to(
            pow  >= 1 - beta,
            toer <= alpha
        ),
        initial_design = OneStageDesign(100, 1.97)
    )

    expect_equal(
        opt_os$design@c1f,
        qnorm(1 - alpha),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        opt_os$design@n1,
        ((qnorm(1 - beta) + qnorm(1 - alpha)) / 0.4)^2,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

}) # end 'optimal one-stage design can be computed'



test_that("Optimal group-sequential design is computable", {

    initial_design_gs <- GroupSequentialDesign(25, 0, 2, 35, 1.96, order)

    opt_gs <<- minimize(
        ess,
        subject_to(
            pow  >= 1 - beta,
            toer <= alpha
        ),
        initial_design_gs
    )

    expect_equal(
        round(evaluate(pow, opt_gs$design), 1),
        0.8,
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        round(evaluate(toer, opt_gs$design), 2),
        0.05,
        tolerance = sqrt(.Machine$double.eps), scale = 1)


    # Check if n2 is equal at boundaries
    expect_equal(
        n2(opt_gs$design, opt_gs$design@c1f),
        n2(opt_gs$design, opt_gs$design@c1e))

    expect_equal(
        opt_gs$nloptr_return$solution[1],
        opt_gs$design@n1)

}) # end 'optimal group-sequential design is computable'



test_that("Optimal group-sequential design is superior to standard gs design", {

    # Create design from rpact
    design_rp <- rpact::getDesignInverseNormal(
        kMax = 2,
        alpha = alpha,
        beta = beta,
        futilityBounds = 0,
        typeOfDesign = "P")
    res <- rpact::getSampleSizeMeans(
        design_rp, normalApproximation = TRUE, alternative = .4 * sqrt(2))
    c2_fun <- function(z){
        w1 <- 1 / sqrt(2)
        w2 <- sqrt(1 - w1^2)
        out <- (design_rp$criticalValues[2] - w1 * z) / w2
        return(out)
    }
    c1f <- qnorm(
        rpact::getDesignCharacteristics(design_rp)$futilityProbabilities
    ) + sqrt(res$numberOfSubjects1[1]) * .4
    rpact_design <- GroupSequentialDesign(
            ceiling(res$numberOfSubjects1[1,]),
        c1f,
        design_rp$criticalValues[1],
        ceiling(res$numberOfSubjects1[2,]),
        rep(2.0, 100),
        100L
    )
    rpact_design@c2_pivots <- sapply(adoptr:::scaled_integration_pivots(rpact_design), c2_fun)

    # use opt_gs from above
    testthat::expect_lte(
        evaluate(ess, opt_gs$design),
        evaluate(ess, rpact_design))

})



test_that("base-case satisfies constraints", {

    opt_ts <<- minimize(
        ess,
        subject_to(
            pow  >= 1 - beta,
            toer <= alpha
        ),
        initial_design
    )

    # compute summaries
    out <- summary(opt_ts$design, "power" = pow, "toer" = toer)

    expect_equal(
        as.numeric(out$scores["power"]),
        0.8,
        tolerance = 1e-3, scale = 1)

    expect_equal(
        as.numeric(out$scores["toer"]),
        0.05,
        tolerance = 1e-3, scale = 1)

}) # end base-case respects constraints



test_that("base-case results are consistent - no post processing", {

    # optimal two-stage design better than optimal group-sequential design
    expect_lt(
        evaluate(ess, opt_ts$design),
        evaluate(ess, opt_gs$design))

    # optimal group-sequential design better than optimal one-stage design
    expect_lt(
        evaluate(ess, opt_gs$design),
        evaluate(ess, opt_os$design))


    # simulate on boundary of null
    sim_null <- simulate(
        opt_ts$design, nsim = 10^6, dist = datadist, theta = .0, seed = 54)

    # check type one error rate on boundary of null
    expect_equal(
        mean(sim_null$reject),
        alpha,
        tolerance = 1e-3, scale = 1)

    # expected sample size on boundary of null
    expect_equal(
        mean(sim_null$n2 + sim_null$n1),
        evaluate(ess_0, opt_ts$design),
        tolerance = 1e-2, scale = 1)

    # simulate under alternative
    sim_alt  <- simulate(
        opt_ts$design, nsim = 10^6, dist = datadist, theta = .4, seed = 54)

    # check power constraint
    expect_equal(
        mean(sim_alt$reject),
        1 - beta,
        tolerance = 1e-2, scale = 1)

    # check expected sample size under alternative
    expect_equal(
        mean(sim_alt$n2 + sim_alt$n1),
        evaluate(ess, opt_ts$design),
        tolerance = .5, scale = 1)

}) # end 'base-case results are consistent'




test_that("conditional constraints work", {

    opt_ts <- suppressWarnings(

        minimize( # ignore: initial design is infeasible
            ess,
            subject_to(
                pow  >= 1 - beta,
                toer <= alpha,
                cp   >= 0.75,
                cp   <= 0.95
            ),
            initial_design,
            opts = list(
                algorithm   = "NLOPT_LN_COBYLA",
                xtol_rel    = 1e-4,
                maxeval     = 10000
            )
        )
    )

    tol <- .005
    # check lower boundary on conditional power
    expect_gte(
        evaluate(cp, opt_ts$design, opt_ts$design@c1f),
        .75 - tol)

    # check lower boundary on conditional power
    expect_lte(
        evaluate(cp, opt_ts$design, opt_ts$design@c1e),
        .95 + tol)

    # test that c2 is monotonously increasing
    expect_true(all(
        sign(diff(opt_ts$design@c2_pivots)) == -1))

}) # end 'conditional constraints work'



test_that("conditional constraints work", {

    expect_equal(
        capture.output(print(opt_os)),
        "optimized OneStageDesign<39;1.64> "
    )


}) # end 'conditional constraints work'
