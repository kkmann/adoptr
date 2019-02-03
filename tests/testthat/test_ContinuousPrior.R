context("ContinuousPrior                                                      ")

test_that("Constructor works", {

    support <- c(0, 1)
    unif <<- ContinuousPrior(function(x) 2*x, support) # define for later use
    expect_equal(
        support,
        unif@support
    )
    expect_equal(
        stats::integrate(unif@pdf, support[1], support[2])$value,
        1
    )

    expect_error(
        # pdf does not integrate to 1
        ContinuousPrior(function(x) x, support)
    )

}) # end 'constructor works'



test_that("bounds() works", {

    expect_equal(
        c(0, 1),
        bounds(unif)
    )

}) # end 'bounds() works'



test_that("expectation() works", {

    expect_equal(
        expectation(unif, identity),
        2/3
    )

}) # end 'expectation works'



test_that("predictive_pdf integrates to 1", {

    normal <<- Normal() # define for later reuse
    n1     <<- 20
    expect_equal(
        stats::integrate(
            function(x1) predictive_pdf(normal, unif, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        1,
        tolerance = .005
    )

}) # end 'predictive_pdf integrates to 1'



test_that("predictive expectation under uniform is larger than 0", {

    expect_gt(
        stats::integrate(
            function(x1) x1 * predictive_pdf(normal, unif, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        0
    )

}) # end 'predictive expectation under uniform is larger than 0'



test_that("conditioning on c(0, .5) leads to correct bounds", {

    unif_cond <<- condition(unif, c(.0, .5))
    expect_equal(
        c(0, .5),
        bounds(unif_cond)
    )

}) # end 'conditioning on c(0, .5) leads to correct bounds'



test_that("conditional predictive pdf integrates to 1", {

    expect_equal(
        stats::integrate(
            function(x1) predictive_pdf(normal, unif_cond, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        1,
        tolerance = .005
    )

}) # end 'conditional predictive pdf integrates to 1'



test_that("conditional unif on c(0, .5) has lower expected value than unconditional", {

    expect_gt(
        stats::integrate(
            function(x1) x1 * predictive_pdf(normal, unif, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        stats::integrate(
            function(x1) x1 * predictive_pdf(normal, unif_cond, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value
    )

}) # end 'conditional unif on c(0, .5) has lower expected value than unconditional'



test_that("posterior pdf integrates to 1", {

    delta <- .5
    x1    <- delta * sqrt(n1)
    post  <<- posterior(normal, unif, x1, n1)

    expect_equal(
        stats::integrate(
            function(theta) post@pdf(theta),
            bounds(post)[1], bounds(post)[2],
            abs.tol = .0001
        )$value,
        1,
        tolerance = .005
    )

}) # end 'posterior pdf integrates to 1'



test_that("bounds of posterior are correct", {

    expect_equal(
        bounds(post),
        bounds(unif)
    )

})



test_that("observing positive z value in normal model results in larger expected value for posterior", {

    expect_gt(
        expectation(post, identity),
        expectation(unif, identity)
    )

}) # end 'observing positive z value in normal model results in larger expected value for posterior'



test_that("increased n lets posterior expectation converge", {

    n1 <- c(10, 20, 33)
    x1 <- .5 * sqrt(n1)

    posteriors <- list()
    for (i in 1:length(n1)) {
        posteriors[[i]] <- posterior(normal, unif, n1[i], x1[i])
    }

    post_expectations <- sapply(posteriors, function(x) expectation(x, identity))

    expect_all( # sequence of posterior expectation should converge to true theta .5
        diff(abs(post_expectations - .5)) < 0
    )

    # TODO: must also work for up to n1 = 100
    # TODO: add test case for two_arm = FALSE!

}) # end '"increased n lets posterior expectation converge"'
