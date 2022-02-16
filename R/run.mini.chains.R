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
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

run.mini.chains <- function(mini.chains, replicates, parallel, path, file.name = "mini.chains", record.tree = FALSE) {

    if(!missing(parallel)) {
        warning("parallel not implemented yet.")
    }

    ## mini.chains
    output <- replicate(replicates, run.a.chain(mini.chains, record.tree), simplify = FALSE)
    class(output) <- c("mini.chains")

    ## Saving the business
    if(!missing(path)) {
        save(output, file = paste0(c(path, file.name, "_", replicates, ".rda"), collapse = "")) #, "_", format(Sys.time(), "%Y-%m-%d-%H%M%S")
    }
    return(output)
}

run.a.chain <- function(mini.chains, record.tree = FALSE) {
    if(!record.tree) {
        return(mini.chains[[sample(1, 1:length(mini.chains))]]$run())
    } else {
        selected_tree <- sample(1, 1:length(mini.chains))
        ## Run the chain
        chain <- mini.chains[[selected_tree]]$run()
        ## Add the tree
        chain$tree <- mini.chains[[selected_tree]]$tree
        return(chain)
    }
}


