## Test
load("MCMCM_chains_test.rda")
chains <- MCMCM_chains_test

test_that("get.burnin works", {
    test <- get.burnin(chains[[1]], buffer = 0)
    expect_equal(test, 700)
    test <- get.burnin(chains[[1]], buffer = 0.25)
    expect_equal(test, 900)
    test <- get.burnin(chains[[2]], buffer = 0)
    expect_equal(test, 400)
    test <- get.burnin(chains[[3]], buffer = 0)
    expect_equal(test, 600)
})

test_that("get.prior works", {

    test <- get.prior(chains[[1]])
    expect_equal(names(test), c("R", "G"))
    expect_equal(names(test[[1]]), c("R1"))
    expect_equal(names(test[[2]]), c("G1", "G2", "G3", "G4"))
    expect_equal(names(test[[1]]$R1), c("V", "nu"))
    expect_equal(names(test[[2]]$G1), c("V", "nu"))
    expect_equal(dim(test[[1]]$R1$V), c(4,4))
    expect_equal(test[[1]]$R1$nu, 0.05)
    expect_equal(round(mean(test[[1]]$R1$V), 6), 0.001468)
    expect_equal(round(mean(test[[2]]$G1$V), 6), 0.029164)
    test <- get.prior(chains[[2]])
    expect_equal(round(mean(test[[1]]$R1$V), 6), 0.001529)
    expect_equal(round(mean(test[[2]]$G1$V), 6), 0.030054)
})

test_that("extract.parameters works", {
    test <- extract.parameters(chains, buffer = 0.25)
    expect_is(test, "list")
    expect_equal(names(test), c("burnin", "priors"))
    expect_equal(test$burnin, 900)
    expect_equal(sum(round(test$priors$G$G1$V, 1)), 0.3)
    expect_equal(sum(round(test$priors$G$G2$V, 1)), 0.2)
    expect_equal(sum(round(test$priors$G$G3$V, 1)), 0.2)
    expect_equal(sum(round(test$priors$G$G4$V, 1)), 0)
    expect_equal(sum(round(test$priors$R$R1$V, 2)), 0.02)
})
