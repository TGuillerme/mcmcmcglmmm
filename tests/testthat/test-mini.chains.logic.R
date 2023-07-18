## Testing the mini chains pipeline
data(morphdat)
data(tree)

test_that("flat.prior works", {

    ## Making a default flat priors for 2 traits
    test <- flat.prior(ntraits = 2)
    expect_is(test, "list")
    expect_equal(names(test), c("R", "G"))
    expect_equal(names(test[[1]]), c("R1"))
    expect_equal(names(test[[2]]), c("G1"))
    expect_equal(dim(test$G$G1$V), c(2,2))
    expect_equal(test$G$G1$nu, 0)

    test <- flat.prior(residuals = 1, randoms = 4, ntraits = 3, nu = 0.002)
    expect_is(test, "list")
    expect_equal(names(test), c("R", "G"))
    expect_equal(names(test[[1]]), c("R1"))
    expect_equal(names(test[[2]]), c("G1", "G2", "G3", "G4"))
    expect_equal(dim(test$G$G1$V), c(3,3))
    expect_equal(test$G$G1$nu, 0.002)
})

test_that("group.terms works", {

    ## Testing clade terms is correct
    groups <- paste("bib", 1:7, sep = ",")
    test <- group.terms(groups, type = "bob")
    expect_is(test, "formula")
    ## Fucking nestedness!
    expect_equal(as.character(test[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[3]]), "bob")
    expect_equal(test[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[3]], 1)
    expect_equal(as.character(test[[2]][[2]][[2]][[2]][[2]][[2]][[3]][[3]]), "bob")
    expect_equal(test[[2]][[2]][[2]][[2]][[2]][[2]][[3]][[2]][[2]][[2]][[3]], 2)
    expect_equal(as.character(test[[2]][[2]][[2]][[2]][[2]][[3]][[3]]), "bob")
    expect_equal(test[[2]][[2]][[2]][[2]][[2]][[3]][[2]][[2]][[2]][[3]], 3)
    expect_equal(as.character(test[[2]][[2]][[2]][[2]][[3]][[3]]), "bob")
    expect_equal(test[[2]][[2]][[2]][[2]][[3]][[2]][[2]][[2]][[3]], 4)
    expect_equal(as.character(test[[2]][[2]][[2]][[3]][[3]]), "bob")
    expect_equal(test[[2]][[2]][[2]][[3]][[2]][[2]][[2]][[3]], 5)
    expect_equal(as.character(test[[2]][[2]][[3]][[3]]), "bob")
    expect_equal(test[[2]][[2]][[3]][[2]][[2]][[2]][[3]], 6)
    expect_equal(as.character(test[[2]][[3]][[3]]), "bob")
    expect_equal(test[[2]][[3]][[2]][[2]][[2]][[3]], 7)
})

test_that("make.mini.chains works", {

    tree_list <- list(tree, tree, tree)
    class(tree_list) <- "multiPhylo"

    ## Model 1.2
    test <- make.mini.chains(data = morphdat, tree = tree_list, dimensions = c(1,2), verbose = FALSE)
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2)trait - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(trait):units")

    ## Model 1
    test <- make.mini.chains(data = morphdat, tree = tree_list, dimensions = c(1,2), verbose = FALSE, residuals = "global")
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2)trait - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(trait):units")


    ## Model 3.2
    test <- make.mini.chains(data = morphdat, tree = tree_list, dimensions = c(1,2), verbose = FALSE, random = "global")
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2)trait - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "~us(trait):animal")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(trait):units")

    ## Model 3
    test <- make.mini.chains(data = morphdat, tree = tree_list, dimensions = c(1,2), verbose = FALSE, random = "global", residuals = "global")
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2)trait - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "~us(trait):animal")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(trait):units")

    ## Model 4
    priors_list <- flat.prior(ntraits = 3, randoms = 1, residuals = 3, nu = 0.1)
    test <- make.mini.chains(data = morphdat, tree = tree, dimensions = c(1:3), verbose = FALSE, residuals = "clade", priors = priors_list, randoms = "global")
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2, PC3)trait:clade - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "~us(trait):animal")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):units + us(at.level(clade, 2):trait):units + us(at.level(clade, 3):trait):units")
    
    ## Model 5
    test <- make.mini.chains(data = morphdat, tree = tree, dimensions = c(1:3), verbose = FALSE, residuals = "clade", randoms = "clade")
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2, PC3)trait:clade - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):animal + us(at.level(clade, 2):trait):animal + us(at.level(clade, 3):trait):animal")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):units + us(at.level(clade, 2):trait):units + us(at.level(clade, 3):trait):units")

    ## Model 6
    test <- make.mini.chains(data = morphdat, tree = tree, dimensions = c(1:3), verbose = FALSE, residuals = "clade", randoms = c("global", "clade"))
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2, PC3)trait:clade - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):animal + us(at.level(clade, 2):trait):animal + us(at.level(clade, 3):trait):animal + us(trait):animal")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):units + us(at.level(clade, 2):trait):units + us(at.level(clade, 3):trait):units")


    ## Model 7
    test <- make.mini.chains(data = morphdat, tree = tree, dimensions = c(1:3), verbose = FALSE, residuals = "global", randoms = c("global", "clade"))
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    tust <- test$run(one_data = test$data, one_tree = test$tree[[1]], params = test$params)
    expect_is(tust, "MCMCglmm")
    expect_equal(paste0(as.character(tust$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2, PC3)trait:clade - 1")
    expect_equal(paste0(as.character(tust$Random$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):animal + us(at.level(clade, 2):trait):animal + us(at.level(clade, 3):trait):animal + us(trait):animal")
    expect_equal(paste0(as.character(tust$Residual$formula), collapse = ""),
                "~us(trait):units")
})

test_that("run/combine.mini.chains works", {

    tree_list <- list(tree, tree, tree)
    class(tree_list) <- "multiPhylo"

    ## Model 3
    mini.chains <- make.mini.chains(data = morphdat, tree = tree_list, dimensions = c(1,2), verbose = FALSE, randoms = "global", residuals = "global")
    expect_is(mini.chains, c("mini.chains"))
    expect_equal(length(mini.chains), 4)
    expect_equal(names(mini.chains), c("data", "tree", "params", "run"))

    ## Works with record tree
    test <- run.mini.chains(mini.chains, replicates = 2, record.tree = TRUE)
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 2)
    expect_is(test[[1]], "MCMCglmm")
    expect_equal(names(test[[1]]), c("Sol","Lambda","ThetaS","VCV","CP","Liab","Fixed","Random","Residual","Deviance","DIC","X","Z","ZR","XL","ginverse","error.term","family","Tune","meta","y.additional","Wscale","tree"))

    ## Test the runnings
    test <- run.mini.chains(mini.chains, replicates = 10)
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 10)
    expect_is(test[[1]], "MCMCglmm")
    expect_equal(names(test[[1]]),c("Sol","Lambda","ThetaS","VCV","CP","Liab","Fixed","Random","Residual","Deviance","DIC","X","Z","ZR","XL","ginverse","error.term","family","Tune","meta","y.additional","Wscale"))
    
    ## Test the combining
    tust <- combine.mini.chains(test)
    expect_is(tust, "MCMCglmm")
    ## Has the correct dimensions for the things that changed
    expect_equal(attr(tust$Sol, "mcpar"), c(11, 811, 10))
    expect_equal(dim(tust$Sol), c(90, 2))
    expect_equal(dim(tust$VCV), c(90, 8))

    ## Test the running
    test2 <- run.mini.chains(mini.chains, replicates = 2, path = "../", file.name = "test_name")

    ## Test the combining
    tust <- combine.mini.chains("../test_name_2.rda")
    expect_is(tust, "MCMCglmm")
    ## Has the correct dimensions for the things that changed
    expect_equal(attr(tust$Sol, "mcpar"), c(11, 171, 10))
    expect_equal(dim(tust$Sol), c(18, 2))
    expect_equal(dim(tust$VCV), c(18, 8))
    expect_true(file.remove("../test_name_2.rda"))
})

test_that("randomised.factors work", {
    data(morphdat)
    data(tree)
    tree_list <- list(tree)
    class(tree_list) <- "multiPhylo"

    priors_list <- flat.prior(ntraits = 3, randoms = 1, residuals = 3, nu = 0.1)
    test <- make.mini.chains(data = morphdat, tree = tree, dimensions = c(1,2,3), verbose = FALSE, residuals = "clade", priors = priors_list, randoms = "global")
    expect_is(test, c("mini.chains"))
    expect_equal(length(test), 4)
    expect_equal(names(test), c("data", "tree", "params", "run"))
    ## Run!
    set.seed(1)
    tust <- run.mini.chains(test, replicates = 1)
    expect_is(tust[[1]], "MCMCglmm")
    expect_equal(paste0(as.character(tust[[1]]$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2, PC3)trait:clade - 1")
    expect_equal(paste0(as.character(tust[[1]]$Random$formula), collapse = ""),
                "~us(trait):animal")
    expect_equal(paste0(as.character(tust[[1]]$Residual$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):units + us(at.level(clade, 2):trait):units + us(at.level(clade, 3):trait):units")
    expect_equal(round(sum(tust[[1]]$Sol), 6), round(5.847793, 6)) 

    ## Same but randomised
    set.seed(1)
    tust <- run.mini.chains(test, replicates = 1, randomised.factors = "clade")
    expect_is(tust[[1]], "MCMCglmm")
    expect_equal(paste0(as.character(tust[[1]]$Fixed$formula), collapse = ""),
                "~cbind(PC1, PC2, PC3)trait:clade - 1")
    expect_equal(paste0(as.character(tust[[1]]$Random$formula), collapse = ""),
                "~us(trait):animal")
    expect_equal(paste0(as.character(tust[[1]]$Residual$formula), collapse = ""),
                "~us(at.level(clade, 1):trait):units + us(at.level(clade, 2):trait):units + us(at.level(clade, 3):trait):units")
    expect_equal(round(sum(tust[[1]]$Sol), 5), round(5.383184, 5))  
})