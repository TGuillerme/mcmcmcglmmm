#' @title run.mini.chains
#'
#' @description Runs a list of mini.chains objects (MCMCglmm)
#'
#' @param mini.chains a \code{"beer"} and \code{"mini.chains"} object.
#' @param replicates the number of replicates per mini chains.
#' @param parallel the number of cores for the paralellisation.
#' @param path optional, the path for saving the data.
#' @param file the prefix for the file name (will be \code{prefix_replicates.rda}).
#' @param record.tree optional, whether to record which tree was used for each replicate (default is \code{FALSE}).
#' @param randomised.factors optional, the names of factors to randomised to create some null models (see details).
#' 
#' @details
#' When using \code{randomise.factors}, you can provide the name of a column in the data (or a vector of columns) to randomise the factors. The factors are randomised so that the same number of elements are contained in each level.
#' If set to \code{NULL} (default), the factors are not randomised.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

run.mini.chains <- function(mini.chains, replicates, parallel, path, file.name = "mini.chains", record.tree = FALSE, randomised.factors = NULL) {

    if(!missing(parallel)) {
        warning("parallel not implemented yet.")
    }

    if(!is.null(randomised.factors)) {
        if(!all(corrects <- randomised.factors %in% colnames(mini.chains$data))) {
            stop(paste0("Factor column named: ", paste(randomised.factors[!corrects], collapse = ", "), " not found in the data."))
        }
    }

    ## Run the mini.chains
    if(!record.tree) {
        output <- replicate(replicates, mini.chains$run(one_data = shuffle.factor(mini.chains$data, randomised.factors), one_tree = mini.chains$tree[[sample(1:length(mini.chains$tree), 1)]], params = mini.chains$params), simplify = FALSE)
    } else {
        ## Variable for saving the recorded tree (out of the environment, hence the long specific name)
        mcmcmcglmmm_recorded_tree_variable <- integer()
        ## Run the mini chains
        output <- replicate(replicates, mini.chains$run(one_data = shuffle.factor(mini.chains$data, randomised.factors), one_tree = mini.chains$tree[[(mcmcmcglmmm_recorded_tree_variable <<- c(sample(1:length(mini.chains$tree), 1), mcmcmcglmmm_recorded_tree_variable))[1]]], params = mini.chains$params), simplify = FALSE)
        ## Save the trees as part of the MCMCglmm objects
        output <- mapply(function(MCMCglmm, tree) {MCMCglmm$tree <- tree; class(MCMCglmm) <- "MCMCglmm"; return(MCMCglmm)}, MCMCglmm = output, tree = trees[mcmcmcglmmm_recorded_tree_variable], SIMPLIFY = FALSE)
    }
    class(output) <- c("mini.chains")

    ## Saving the business
    if(!missing(path)) {
        save(output, file = paste0(c(path, file.name, "_", replicates, ".rda"), collapse = "")) #, "_", format(Sys.time(), "%Y-%m-%d-%H%M%S")
    }
    return(output)
}

shuffle.factor <- function(data, randomised.factors) {
    if(is.null(randomised.factors)) {
        return(data)
    } else {
        for(one_column in 1:length(randomised.factors)) {
            ## Resample the factor
            data[, randomised.factors[one_column]] <- sample(data[, randomised.factors[one_column]], size = nrow(data), replace = FALSE)
        }
        return(data)
    }
}

