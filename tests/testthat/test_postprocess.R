context("check postprocess()")

n1     <- 25
c1f    <-   .0
c1e    <-  2.0
order  <- 5L
n2_piv <- 40.0
c2_piv <- 1.96

initial_design <- TwoStageDesign(n1, c1f, c1e, n2_piv, c2_piv, order)

dist        <- Normal()
null        <- PointMassPrior(.0, 1)
alternative <- PointMassPrior(.4, 1)

ess  <- expected(ConditionalSampleSize(dist, alternative))
pow  <- expected(ConditionalPower(dist, alternative))
toer <- expected(ConditionalPower(dist, null))


test_that("post-processing yields integer sample sizes", {

    res <- suppressWarnings(minimize( # we do not warning about non-convergence!

        ess,

        subject_to(
            pow  >= 0.8,
            toer <=  .05
        ),

        initial_design = initial_design,

        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-4,
            maxeval     = 100 # only the principle is tested
        )

    ))

    res_post <- suppressWarnings(postprocess(res))

    expect_warning(
        postprocess(res)
    )  # warning because maxeval was reached

    # n1 is integer
    expect_equal(
        res_post$design@n1,
        round(res$design@n1)
    )

    # n2 is integer
    expect_equal(
        res_post$design@n2_pivots,
        round(res$design@n2_pivots)
    )

}) # end 'post-processing yields integer sample sizes'



test_that("post-processing of OneStageDesign yields integer sample sizes", {

    res <- suppressWarnings(minimize( # we do not warning about non-convergence!

        ess,

        subject_to(
            pow  >= 0.8,
            toer <=  .05
        ),

        initial_design        = OneStageDesign(50.0, 2.0),

        opts = list(
            algorithm   = "NLOPT_LN_COBYLA",
            xtol_rel    = 1e-4,
            maxeval     = 10 # only the principle is tested
        )

    ))

    res_post <- suppressWarnings(postprocess(res))

    expect_warning(
        postprocess(res)
    )  # warning because maxeval was reached


    # n1 is integer
    expect_equal(
        res_post$design@n1,
        round(res$design@n1)
    )


}) # end 'post-processing of OneStageDesign yields integer sample sizes'
