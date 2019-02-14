context("TwoStageDesign")

test_that("Optimal design with point prior is computable", {
    # define an initial design
    n1     <- 150
    c1f    <-  0.7
    c1e    <-  2.5
    number_knots <- 5L
    n2_piv <- rep(150.0, number_knots)
    c2_piv <- rep(1.96, number_knots)
    design <- gq_design(n1, c1f, c1e, n2_piv, c2_piv, number_knots)

    # check if functions are defined correctly
    expect_equal(
        n2(design, 1.0),
        150.0
    )

    expect_equal(
        c2(design, 1.0),
        1.96
    )

    # check if length does fit
    expect_equal(
        length(tunable_parameters(design)),
        2 * number_knots + 3
    )

    # check if key figures can be computed
    null        <- PointMassPrior(.0, 1)
    alternative <- PointMassPrior(.3, 1)

    dist <- Normal(two_armed = T)

    ess  <- integrate(ConditionalSampleSize(dist, alternative))
    cp   <- ConditionalPower(dist, alternative)
    pow  <- integrate(cp)
    toer <- integrate(ConditionalPower(dist, null))
    smth <- AverageN2()

    expect_equal(
        round(evaluate(ess, design), 1),
        214.8
    )

    expect_equal(
        round(evaluate(pow, design), 3),
        0.858
    )

    expect_equal(
        round(evaluate(toer, design), 3),
        0.012
    )

    #compute optimal design

    objective <- function(x) {
        d  <- update(design, x)
        evaluate(ess, d) + .0001*evaluate(smth, d)
    }

    constraint <- function(x) {
        d  <- update(design, x)
        c(
            .8 - evaluate(pow, d),
            evaluate(toer, d) - 0.025,
            x[2] - x[3] + .1,
            diff(c2(d, scaled_integration_pivots(d)))
        )
    }

    ub <- c(200, 1, 4, numeric(number_knots) + 200, numeric(number_knots) + 3)
    lb <- c(10, -1, 1, numeric(number_knots) + 2, numeric(number_knots) - 3)

    res <- nloptr::nloptr(
        x0 = tunable_parameters(design),
        lb = lb,
        ub = ub,
        eval_f      = objective,
        eval_g_ineq = constraint,
        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-4,
            maxeval     = 5000
        )
    )

    d2 <<- update(design, res$solution)

    expect_equal(
        round(evaluate(pow, d2), 2),
        0.80
    )

    expect_equal(
        round(evaluate(toer, d2), 3),
        0.025
    )

    expect_equal(
        sign(evaluate(ess, d2) - evaluate(ess, design)),
        -1
    )

    expect_equal(
        round(d2@n1),
        100
    )

    expect_equal(
        round(d2@c1e, 1),
        2.3
    )

    expect_equal(
        round(d2@c1f, 1),
        0.8
    )


}) # end 'Optimal design with point prior is computable'


test_that("Optimal design is superior to standard GS design", {

    # Create design from rpact
    design_rp <- rpact::getDesignInverseNormal(
        kMax = 2,
        alpha = 0.025,
        beta = 0.2,
        futilityBounds = 0,
        typeOfDesign = "P"
    )

    res <- rpact::getSampleSizeMeans(
        design_rp, normalApproximation = TRUE, alternative = .3
    )

    char <- rpact::getDesignCharacteristics(design_rp)

    n1 <- res$numberOfPatientsGroup1[1,]
    n2 <- res$numberOfPatientsGroup1[2,]


    f <- function(z){
        w1 <- 1 / sqrt(2)
        w2 <- sqrt(1 - w1^2)
        out <- (design_rp$criticalValues[2] - w1 * z) / w2
        return(out)
    }

    c1f <- qnorm(char$futilityProbabilities) +
        sqrt(res$numberOfPatientsGroup1[1]) * .3 / sqrt(2)
    c1e <- design_rp$criticalValues[1]

    x <- GaussLegendreRule(5)$nodes
    h <- (c1e - c1f) / 2
    x <- h * x + (h + c1f)

    design_gs <- gq_design(
        ceiling(n1),
        c1f,
        c1e,
        rep(ceiling(n2), 5),
        sapply(seq(c1f, c1e, length.out = 5), f),
        5L
    )

    # Define key figures
    ess   <- integrate(ConditionalSampleSize(Normal(), PointMassPrior(.3, 1)))
    pow   <- integrate(ConditionalPower(Normal(), PointMassPrior(.3, 1)))
    toer  <- integrate(ConditionalPower(Normal(), PointMassPrior(.0, 1)))


    expect_gt(
        evaluate(ess, design_gs),
        evaluate(ess, d2)
    )

    expect_equal(
        evaluate(pow, d2),
        evaluate(pow, design_gs),
        tolerance = .01
    )

    expect_equal(
        evaluate(toer, d2),
        evaluate(toer, design_gs),
        tolerance = .005
    )

}) # end 'Optimal design is superior to standard GS design'



test_that("errors are returned correctly", {
    expect_error(
        GaussLegendreRule(-1)
    )

    f <- function(x) x

    expect_error(
        integrate_rule(f, 0, 1, .5, c(1, 1))
    ) # pivots and weigths of same length

    expect_error(
        integrate_rule(f, 0, 3, c(1, 2), c(1, 1))
    ) # x is scaled automatically

    expect_error(
        integrate_rule(f, 0, 1, c(.3, .6), c(0, 1))
    ) # weights must be positive

    order = 3L
    expect_error(
        gq_design(50, 0, 2, rep(50, 2), rep(3, order), order)
    ) # parameters length must fit

}) # end 'errors are returned correctly'




test_that("print methods", {
    vdiffr::expect_doppelganger(
        "Design print",
        print.TwoStageDesignSummary(summary(d2))
    )

})

