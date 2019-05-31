context("ContinuousPrior")



test_that("Constructor works", {

    support <- c(0, 1)
    prior   <<- ContinuousPrior(function(x) 2*x, support) # define for later use

    expect_equal(
        support,
        prior@support,
        tolerance = 1e-6, scale = 1)

    expect_equal(
        stats::integrate(prior@pdf, support[1], support[2])$value,
        1,
        tolerance = 1e-6, scale = 1)

    expect_error(
        # pdf does not integrate to 1
        ContinuousPrior(function(x) x, support),
        tolerance = 1e-6, scale = 1)

}) # end 'constructor works'



test_that("bounds() works", {

    expect_equal(
        c(0, 1),
        bounds(prior),
        tolerance = 1e-6, scale = 1)

}) # end 'bounds() works'



test_that("expectation() works", {

    expect_equal(
        expectation(prior, identity),
        2/3,
        tolerance = 1e-6, scale = 1)

}) # end 'expectation works'



test_that("predictive_pdf integrates to 1", {

    normal <<- Normal() # define for later reuse
    n1     <<- 20
    expect_equal(
        stats::integrate(
            function(x1) predictive_pdf(normal, prior, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        1,
        tolerance = 1e-4, scale = 1)

}) # end 'predictive_pdf integrates to 1'



test_that("predictice_cdf is monotonously increasing", {

    pcdf <- function(x1) predictive_cdf(normal, prior, x1, n1)
    x <- seq(0.0, 2.0, length.out = 10)
    y <- sapply(x, pcdf)

    expect_equal(
        sign(diff(y)),
        rep(1, length(y) - 1),
        tolerance = 1e-6, scale = 1)

}) # end 'predicitive_cdf is monotonously increasing'



test_that("predictive expectation under prior is larger than 0", {

    expect_gt(
        stats::integrate(
            function(x1) x1 * predictive_pdf(normal, prior, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        0)

}) # end 'predictive expectation under prior is larger than 0'



test_that("conditioning works", {

    prior_cond <<- condition(prior, c(.0, .5))
    expect_equal(
        c(.0, .5),
        bounds(prior_cond)
    ) # conditioning on c(0, .5) leads to correct bounds

    unif_prior <- ContinuousPrior(function(x) dunif(x, 0, 1), c(0, 1))
    cond_unif  <- condition(unif_prior, c(-10, 10))
    expect_equal(
        bounds(unif_prior),
        bounds(cond_unif)
    ) # conditioning on inverval larger than the support does not change anything

    expect_length(
        cond_unif@pdf(c(.5, 2)),
        2
    ) # conditioned pdf is vectorized

    expect_error(
        condition(unif_prior, c(2, 3))
    )

    expect_equal(
        cond_unif@pdf(c(-.5, 1.5)),
        c(0, 0)
    ) # pdf is 0 outside support after conditioning

}) # end 'conditioning works'



test_that("tightening can decrease support", {

    support <- c(-5, 5)
    prior_2 <- ContinuousPrior(function(x) dnorm(x, sd = .1), support, tighten_support = TRUE)

    expect_gt(
        prior_2@support[1],
        support[1])

    expect_lt(
        prior_2@support[2],
        support[2])

}) # end 'tightening can decrease support'



test_that("conditional predictive pdf integrates to 1", {

    expect_equal(
        stats::integrate(
            function(x1) predictive_pdf(normal, prior_cond, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        1,
        tolerance = 1e-4, scale = 1)

}) # end 'conditional predictive pdf integrates to 1'



test_that("conditional prior on c(0, .5) has lower expected value than unconditional", {

    expect_gt(
        stats::integrate(
            function(x1) x1 * predictive_pdf(normal, prior, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value,
        stats::integrate(
            function(x1) x1 * predictive_pdf(normal, prior_cond, x1, n1),
            qnorm(.0005), qnorm(.9995, mean = sqrt(n1)),
            abs.tol = .0001
        )$value)

}) # end 'conditional prior on c(0, .5) has lower expected value than unconditional'



test_that("posterior pdf integrates to 1", {

    delta <- .5
    x1    <- delta * sqrt(n1)
    post  <<- posterior(normal, prior, x1, n1)

    expect_equal(
        stats::integrate(
            function(theta) post@pdf(theta),
            bounds(post)[1], bounds(post)[2],
            abs.tol = .0001
        )$value,
        1,
        tolerance = 1e-4, scale = 1)

}) # end 'posterior pdf integrates to 1'



test_that("bounds of posterior are correct", {

    expect_equal(
        bounds(post),
        bounds(prior),
        tolerance = 1e-4, scale = 1)

})



test_that("observing positive z value in normal model results in larger expected value for posterior", {

    expect_gt(
        expectation(post, identity),
        expectation(prior, identity))

}) # end 'observing positive z value in normal model results in larger expected value for posterior'



test_that("increased n lets posterior expectation converge", {

    # one-arm case
    n1         <- c(10, 20, 33, 100)
    x1         <- .5 * sqrt(n1)
    normal     <- Normal(two_armed = FALSE)
    posteriors <- list()
    for (i in 1:length(n1)) {
        posteriors[[i]] <- posterior(normal, prior, x1[i], n1[i])
    }
    post_expectations <- sapply(posteriors, function(x) expectation(x, identity))

    expect_true( # sequence of posterior expectation should converge to true theta .5
        all(diff(abs(post_expectations - .5)) < 0))

    # two-arm case
    n1         <- c(10, 20, 33, 250)
    x1         <- .5 * sqrt(n1) / sqrt(2)
    normal     <- Normal(two_armed = TRUE)
    posteriors <- list()
    for (i in 1:length(n1)) {
        posteriors[[i]] <- posterior(normal, prior, x1[i], n1[i])
    }
    post_expectations <- sapply(posteriors, function(x) expectation(x, identity))

    expect_true( # sequence of posterior expectation should converge to true theta .5
        all(diff(abs(post_expectations - .5)) < 0))

}) # end '"increased n lets posterior expectation converge"'



test_that("Errors are defined correctly", {

    expect_error(
        ContinuousPrior(function(x) 2*x, 1)) # support must be of length 2

    expect_error(
        ContinuousPrior(function(x) 1/x^2, c(1, Inf))) # support must be finite

    expect_error(
        ContinuousPrior(function(x) 2*x, c(1, 0))) # support[2] must be greater or equal than support[1]

    # the same in 'condition'
    prior <- ContinuousPrior(function(x) dnorm(x, mean = 0, sd = .1), c(-5, 5))

    expect_error(
        condition(prior, 1)) # interval must be of length 2

    expect_error(
        condition(prior, c(0, Inf))) # interval must be finite

    expect_error(
        condition(prior, c(3, 0))) # interval[2] must be greater or equal than interval[1]

    condprior <- condition(prior, c(0,3))

    expect_equal(
        stats::integrate(function(x) condprior@pdf(x), 0, 3)$value,
        1,
        tolerance = 1e-4, scale = 1) # conditioning works when defined correctly

    expect_error(
        posterior(Normal(), prior, c(1,2), 50)) # posterior not vectorized in x1

}) # end 'errors are defined correctly'



test_that("show method returns class name", {

    prior <- ContinuousPrior(function(x) 2*x, c(0, 1))

    expect_true(
        class(prior)[1] == capture.output(show(prior)))

}) # end 'show method returns class name'
