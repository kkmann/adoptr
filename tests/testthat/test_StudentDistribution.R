context("Student data distribution")

library(adoptr)
library(testthat)

test_that("Constructor works", {

  expect_true(
    Student()@two_armed)

  expect_true(
    !Student(two_armed = FALSE)@two_armed)

  expect_equal(
    Student(TRUE)@multiplier, 2)

  expect_equal(
    Student(FALSE)@multiplier, 1)

}) # end 'constructor works'


test_that("pdf is defined correctly", {

  dist    <- Student(TRUE)
  x       <- seq(-3, 3, by = .1)
  n       <- 22
  theta   <- .35
  expect_equal(
    probability_density_function(dist, x, n, theta),
    stats::dt(x,  df = 2 * (n - 1), ncp = sqrt(n / 2) * theta),
    tolerance = 1e-6, scale = 1)
})



test_that("cdf is defined correctly", {

  dist    <- Student(FALSE)
  x       <- seq(-3, 3, by = .1)
  n       <- 11
  theta   <- 0.24
  expect_equal(
    cumulative_distribution_function(dist, x, n, theta),
    stats::pt(x,  df = 1 * (n - 1), ncp = sqrt(n / 1) * theta),
    tolerance = 1e-6, scale = 1)
})



test_that("quantile is defined correctly", {

  dist  <- Student(TRUE)
  probs <- seq(.01, 1, by = .01)
  n     <- 7
  theta <- -.1
  expect_equal(
    quantile(dist, probs, n, theta),
    stats::qt(probs, df = 2 * (n - 1), ncp = sqrt(n / 2) * theta),
    tolerance = 1e-6, scale = 1)

  expect_warning(
    quantile(dist, -1, n, theta)
  )

})


test_that("vectorization works", {
  dist  <- Student(two_armed = TRUE)
  prior <- PointMassPrior(c(.3, .4), c(.4, .6))
  expect_equal(
    predictive_pdf(dist, prior, c(0, 1), 10),
    sapply(c(0, 1), function(x) predictive_pdf(dist, prior, x, 10))
  )
})


test_that("simulate respects seed", {

  expect_equal(
    simulate(Student(), 10, 10, 0, seed = 42),
    simulate(Student(), 10, 10, 0, seed = 42),
    tolerance = 1e-6, scale = 1)

  set.seed(42)

  expect_true(
    all(simulate(Student(), 10, 12, -.1) != simulate(Student(), 10, 12, -.1)))

}) # end 'simulate respects seed'



test_that("Student needs larger sample sizes than normal", {
  dist1 <- Student()
  pow1  <- Power(dist1, PointMassPrior(.3,1))
  ess1  <- ExpectedSampleSize(dist1, PointMassPrior(.3, 1))
  toer1 <- Power(dist1, PointMassPrior(0,1))
  opt1  <- minimize(ess1, subject_to(pow1>=.9, toer1<=.025), OneStageDesign(200,2))$design

  dist2 <- Normal()
  pow2  <- Power(dist2, PointMassPrior(.3,1))
  ess2  <- ExpectedSampleSize(dist2, PointMassPrior(.3, 1))
  toer2 <- Power(dist2, PointMassPrior(0,1))
  opt2  <- minimize(ess2, subject_to(pow2>=.9, toer2<=.025), OneStageDesign(200,2))$design

  expect_lte(
    opt2@n1,
    opt1@n1
  )
}) # end 'student needs larger sample sizes than normal'


test_that("Published results are reproducible", {
  # source: Chapter 3 of M. Kieser (2020):
  # Methods and Applications of Sample Size Calculation and Recalculation in Clinical Trials

  dist <- Student(TRUE)
  H_0  <- PointMassPrior(0, 1)
  ess  <- ExpectedSampleSize(dist, H_0)
  toer <- Power(dist, H_0)


  pow1  <- Power(dist, PointMassPrior(1.5, 1))
  osd1  <- minimize(ess, subject_to(toer <= 0.025, pow1 >= 0.8),
                    OneStageDesign(20, 2))$design
  expect_equal(ceiling(osd1@n1), 18 / 2)

  pow2  <- Power(dist, PointMassPrior(1.0, 1))
  osd2  <- minimize(ess, subject_to(toer <= 0.025, pow2 >= 0.8),
                    OneStageDesign(20, 2))$design
  expect_equal(ceiling(osd2@n1), 34 / 2)

  pow3  <- Power(dist, PointMassPrior(0.5, 1))
  osd3  <- minimize(ess, subject_to(toer <= 0.025, pow3 >= 0.8),
                    OneStageDesign(20, 2))$design
  expect_equal(ceiling(osd3@n1), 128 / 2)


}) # end 'published results are reproducible'




test_that("show method", {

  expect_equal(
    capture.output(show(Student())),
    "Student<two-armed> "
  )

  expect_equal(
    capture.output(show(Student(two_armed = FALSE))),
    "Student<single-armed> "
  )

})
