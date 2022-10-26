context("check minimize()")



# preliminaries
order <- 5L

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

initial_design <- get_initial_design(.4, alpha, beta, "two-stage", dist=datadist, order=order)



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
        initial_design = get_initial_design(.4, alpha, beta, "one-stage", dist=datadist, order=order)
    )

    expect_equal(
        opt_os$design@c1f,
        qnorm(1 - alpha),
        tolerance = sqrt(.Machine$double.eps), scale = 1)

    expect_equal(
        opt_os$design@n1,
        ((qnorm(1 - beta) + qnorm(1 - alpha)) / 0.4)^2,
        tolerance = 1e-4, scale = 1)

}) # end 'optimal one-stage design can be computed'



test_that("Optimal group-sequential design is computable", {

    initial_design_gs <- get_initial_design(.4, alpha, beta, "group-sequential", dist=datadist, order=order)

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
    out <- summary(opt_ts$design, "power" = pow, "toer" = toer, "CP" = cp, rounded = FALSE)

    out2 <- summary(opt_ts$design, "power" = pow, "toer" = toer, rounded = FALSE)

    out3 <- summary(opt_ts$design,  "CP" = cp, rounded = FALSE)

    out4 <- summary(opt_ts$design, rounded = TRUE)

    expect_equal(
        as.numeric(out$uncond_scores["power"]),
        0.8,
        tolerance = 1e-3, scale = 1)

    expect_equal(
        as.numeric(out$uncond_scores["power"]),
        as.numeric(out2$uncond_scores["power"]),
        tolerance = 1e-3, scale = 1)

    expect_equal(
        as.numeric(out$uncond_scores["toer"]),
        0.05,
        tolerance = 1e-3, scale = 1)


    expect_equal(
        as.numeric(out$uncond_scores["toer"]),
        as.numeric(out2$uncond_scores["toer"]),
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

    # maximum sample size of adaptive design is larger than of one-stage design
    mss <- MaximumSampleSize()
    expect_lte(
        evaluate(mss, opt_os$design),
        evaluate(mss, opt_ts$design)
    )

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
        "OneStageDesign<optimized;n=39;c=1.64> "
    )


}) # end 'conditional constraints work'




test_that("initial design works", {
    expect_error(
        get_initial_design(.4, .025, .2, "adaptive", dist=Normal(), order=6L)
    )

    expect_error(
        get_initial_design(0.4,0.025,0.2, dist=Survival(0.8))
    )

    expect_error(
        get_initial_design(.4, 1.025, .2, "two-stage", dist=Normal(), order=6L)
    )

    expect_true(
        is(get_initial_design(.4, .025, .2, "two-stage", dist=Normal(F), order=6L), "TwoStageDesign")
    )

    expect_true(
        is(get_initial_design(1.4, .025, .2, "two-stage", type_n2="linear_decreasing",dist=Survival(0.7), order=6L,slope=-10),
           "TwoStageDesignSurvival")
    )

    expect_true(
        is(get_initial_design(1.4, .025, .2, "two-stage", dist=Survival(0.7), order=6L),
           "TwoStageDesignSurvival")
    )

    expect_true(
        is(get_initial_design(1.4, .025, .2, "two-stage", type_n2 = "linear_increasing",dist=Survival(0.7), order=6L,slope=10),
           "TwoStageDesignSurvival")
    )

    expect_true(
        is(get_initial_design(.4, .025, .2, "group-sequential", dist=Normal(), order=6L), "GroupSequentialDesign")
    )

    expect_true(
        is(get_initial_design(.4, .025, .2, "group-sequential", dist=Normal(), type_c2 = "constant",
                              order=6L), "GroupSequentialDesign")
    )

    expect_true(
        is(get_initial_design(1.4, .025, .2, "group-sequential", dist=Survival(0.7), order=6L),
           "GroupSequentialDesignSurvival")
    )

    expect_true(
        is(get_initial_design(.4, .025, .2, "one-stage", dist=Normal(), order=6L), "OneStageDesign")
    )


    init <- get_initial_design(.4, .025, .2, "two-stage", dist=Normal(F), cf=0.5,
                          type_c2="constant", ce=2.5)

    #check values futility boundary
    expect_true(
        init@c1f == 0.5
    )

    # check that c2 is constant
    expect_true(
        all(init@c2_pivots-init@c2_pivots[1]==rep(0,7))
    )

    #check efficacy boundary
    expect_true(
        init@c1e==2.5
    )

    init2 <- get_initial_design(.4, .025, .2, "two-stage", dist=Normal(T),
                                type_n2="linear_increasing",slope=20)

    #check that n2 is increasing
    expect_false(
        is.unsorted(init2@n2_pivots)
    )

    init2new <- get_initial_design(.4, .025, .2, "two-stage", dist=Normal(T),
                                type_n2="linear_increasing", info_ratio = 0.3)

    #check that n2 is increasing
    expect_false(
        is.unsorted(init2new@n2_pivots)
    )


    #check the slope
    expect_true(
        n2(init2,1.5)-n2(init2,0.5) == 20
    )

    init3 <- get_initial_design(.4, .025, .2, "two-stage", dist=Normal(T),
                                type_c2="constant", info_ratio = 0.2,ce=2)

    #check information ratio
    expect_equal(4*rep(init3@n1,7),init3@n2_pivots)


    #check that student distribution uses normal distribution
    init4 <- get_initial_design(0.4,0.025,0.1,"one-stage",dist=Student())
    expect_equal(
        init4@c1e,
        1.96,
        tolerance=1e-3
    )

    #check warnings are raised when unnecessarily specifying n2 or c2 type
    expect_warning(
        get_initial_design(0.4,0.025,0.2,"one-stage",type_c2 = "constant")
    )

    expect_warning(
        get_initial_design(0.4,0.025,0.2,"one-stage",type_n2 = "constant")
    )

    #check warning when ce is chosen too small

    expect_warning(
        get_initial_design(0.4,0.025,0.2,"two-stage",type_c2="constant",ce=1)
    )

    expect_warning(
        get_initial_design(0.4,0.025,0.2,"two-stage",type_c2="linear_decreasing",ce=1)
    )

    #check that ce can be chosen

    expect_equal(
        get_initial_design(0.3,0.025,0.2,type_c2="linear_decreasing",ce=3)@c1e,
        3
    )

    init5 <- get_initial_design(.4, .025, .2, "two-stage", dist=Normal(F), cf=0.5, type_n2 = "constant",
                       type_c2="constant", ce=2.5)

    #check that c2 is constant
    expect_true(
            all(init@c2_pivots-init@c2_pivots[1]==rep(0,7))
        )

    #check that c2 is linear decreasing
    init6 <- get_initial_design(1.4, .025, .2, "two-stage", dist=Survival(0.7), cf=0.5, type_n2 = "constant",
                                type_c2="linear_decreasing")

    expect_false(
        is.unsorted(rev(init6@c2_pivots))
    )

    #check linear decreasing n2 design
    init7 <- get_initial_design(0.4, .025, .2, "two-stage", cf=0.5,
                                type_n2="linear_decreasing")

    expect_false(
        is.unsorted(rev(init7@n2_pivots))
    )

    expect_error(
        get_initial_design(0.4, .025, .2, "two-stage", cf=0.5,
                           type_n2="linear_decreasing", slope=10)
    )

    expect_error(
        get_initial_design(0.4, .025, .2, "two-stage", cf=0.5,
                           type_n2="linear_increasing", slope=-10)
    )


    #check that slope is automatically adjusted
    expect_warning(
        get_initial_design(0.4, .025, .2, "two-stage", cf=0.5,
                                 type_n2="linear_decreasing",slope=-100))

    expect_warning(
        get_initial_design(0.4, .025, .2, "two-stage", cf=0.5,
                           type_n2="linear_increasing",slope=100))



}) # end 'initial design works'


test_that("constraint checks are working", {
    init <- TwoStageDesign(25,0.2,2.5,c(80,77,70,61,50,36,25),c(2.2,2.1,1.8,1.4,0.9,0.3,-0.1),order=7L)

    expect_warning(
        opt_design <- minimize(ess,subject_to(toer<=0.01,pow>=0.95,cp>=0.95,MaximumSampleSize()<=70),
                               initial_design = init, check_constraints=TRUE,
                               opts=list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 1))

    )

    expect_warning(
        opt_design <- minimize(ess,subject_to(toer<=0.025,pow>=0.8,MaximumSampleSize()<=pow,cp>=ConditionalSampleSize()),
                               initial_design = init, check_constraints = TRUE,
                               opts=list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 1))
    )
}) # end 'constraint checks are working'
