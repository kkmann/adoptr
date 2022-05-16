context("Survival class")

test_that("Constructor works",  {

    expect_true(
        Survival()@two_armed)

    expect_true(
        !Survival(two_armed=FALSE)@two_armed)

    expect_error(
        Survival(-.1))

    expect_error(
        Survival(1.1))
})

test_that("pdf is defined correctly" , {

    dist <- Survival(0.7)


})
