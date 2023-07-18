#' @title Extract parameters
#'
#' @description Extracting parameters from prior MCMCglmm chains
#'
#' @param chains The prior MCMCglmm chains
#' @param parameters The parameters to extract (see details).
#' @param buffer The buffer for the burnin (how many itteration to include past the first median estimate)
#' @param nu The degree of belief parameter (nu) for the priors.
#' 
#' @details
#' The parameters that can be extracted are:
#' \itemize{
#'      \item \code{"burnin"}: the global burnin length
#'      \item \code{"priors"}: a list of priors (from the posteriors)
#' }
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
extract.parameters <- function(chains, parameters = c("burnin", "priors"), buffer = 0.1, nu = 0.05) {

    ## Get the required parameters
    param_out <- list()
    if("burnin" %in% parameters) {
        param_out$burnin <- max(unlist(lapply(chains, get.burnin, buffer = buffer)))
    }
    if("priors" %in% parameters) {
        ## Check for burnin
        if(!is.null(param_out$burnin)) {
            burnin_n <- param_out$burnin/attr(chains[[1]]$Sol, "mcpar")[3]
        } else {
            burnin_n <- max(unlist(lapply(chains, get.burnin, buffer = buffer)))/attr(chains[[1]]$Sol, "mcpar")[3]
        }
        ## Get the list of priors
        priors_list <- lapply(chains, get.prior, nu = nu, burnin_n = burnin_n)
        ## Merging the priors
        param_out$priors <- merge.V(priors_list)
    }

    return(param_out)
}

## Merge the lists of priors V matrices
merge.V <- function(priors_list) {
    ## Get the number of levels
    n_ran <- length(priors_list[[1]]$G)
    n_res <- length(priors_list[[1]]$R)
    n_traits <- dim(priors_list[[1]]$R$R1$V)[1]
    nu <- priors_list[[1]]$R$R1$nu

    ## Get the template to fill
    template <- flat.prior(ntraits = n_traits, residuals = n_res, randoms = n_ran, nu = nu)
    
    ## Get the length of the list
    n_chains <- length(priors_list)

    ## Fill the list for G
    if(n_ran != 0) {
        for(one_ran in 1:n_ran) {
            ## Get the list of V matrices
            V_list <- lapply(lapply(priors_list, `[[`, "G"), function(X, i) return(X[[i]]$V), i = one_ran)
            ## Get the prior matrix
            prior_matrix <- positive.definite(V_list)
            if(is.null(prior_matrix)) {
                stop(paste0("Impossible to get a positive definite matrix from the posteriors of the random term number ", one_ran, "."), call. = FALSE)
            } else {
                if(is(prior_matrix, "list")) {
                    ## Warning
                    warning(paste0("The median of the posteriors for the random term number ", one_ran, " didn't resulted in a positive definite matrix. The mean was used instead."), call. = FALSE)
                    prior_matrix <- prior_matrix$matrix
                }
            }

            ## Refill the template
            template$G[[one_ran]]$V <- prior_matrix
        }
    }

    ## Fill the list for R
    if(n_res != 0) {
        for(one_res in 1:n_res) {
            ## Get the list of V matrices
            V_list <- lapply(lapply(priors_list, `[[`, "R"), function(X, i) return(X[[i]]$V), i = one_res)

            ## Get the prior matrix
            prior_matrix <- positive.definite(V_list)
            if(is.null(prior_matrix)) {
                stop(paste0("Impossible to get a positive definite matrix from the posteriors of the residual term number ", one_res, "."), call. = FALSE)
            } else {
                if(is(prior_matrix, "list")) {
                    ## Warning
                    warning(paste0("The median of the posteriors for the residual term number ", one_res, " didn't resulted in a positive definite matrix. The mean was used instead."), call. = FALSE)
                    prior_matrix <- prior_matrix$matrix
                }
            }

            ## Refill the template            
            template$R[[one_res]]$V <- apply(array(do.call(cbind, V_list), dim = c(dim(V_list[[1]]), length(V_list))), c(1,2), median)
        }
    }
    return(template)
}

## Make sure the matrix is positive definite
positive.definite <- function(V_list) {

    ## Get the median
    matrix <- apply(array(do.call(cbind, V_list), dim = c(dim(V_list[[1]]), length(V_list))), c(1,2), median)

    ## Check if the matrix is positive definitive (all eigen values > 0)
    if(any((eig_val <- eigen(matrix)$values) < 0)) {
        ## Use the mean rather than the median
        matrix <- apply(array(do.call(cbind, V_list), dim = c(dim(V_list[[1]]), length(V_list))), c(1,2), mean)
        if(any((eig_val <- eigen(matrix)$values) < 0)) {
            matrix <- NULL
        } else {
            matrix <- list(matrix = matrix, use_mean = TRUE)
        }
    }
    return(matrix)
}

## Get the burnin value
get.burnin <- function(chain, buffer) {

    ## Find the median for each MCMC
    find.median <- function(estimate, buffer) {
        ## Get the first estimate past the median
        return(ceiling(which(estimate > median(estimate))[1] * (1+buffer)))
    }

    burnin_points <- apply(chain$Sol, 2, find.median, buffer)
    return(max(burnin_points * attr(chain$Sol, "mcpar")[3]))
}

## Get the priors from one chain
get.prior <- function(chain, nu = 0.05, burnin_n) {

    ## Get the set of parameters
    traits <- MCMCglmm.traits(chain)
    n_traits <- length(traits)
    ## The number of levels
    levels <- MCMCglmm.levels(chain)
    raw_levels <- MCMCglmm.levels(chain, convert = FALSE)
    n_levels <- length(levels)
    n_ran <- sum(names(levels) %in% "random")
    n_res <- sum(names(levels) %in% "residual")

    ## Create the prior template
    template <- flat.prior(ntraits = n_traits, residuals = n_res, randoms = n_ran, nu = nu)

    ## Get the burnin
    if(missing(burnin_n)) {
        burnin_n <- as.integer(get.burnin(chain, buffer = 0.25)/attr(chain$Sol, "mcpar")[3])
    }
    burnin_n <- 1:burnin_n

    ## Handle the G-Structure bits
    if(n_ran != 0) {
        for(i in 1:n_ran) {
            ## Get the columns of interest
            #cells <- grep(gsub(")$", "",gsub(":animal", "", gsub("us\\(", "", raw_levels[(1:n_ran)[i]]))), colnames(chain$VCV), fixed = TRUE)
            cells <- (1:n_traits^2)+(n_traits^2)*(i-1)

            ## Get the mean estimates
            template$G[[i]]$V <- matrix(apply(chain$VCV[-burnin_n, cells], 2, mean), n_traits, n_traits, byrow = FALSE)
        }
    }

    ## Handle the R-Structure bits
    if(n_res != 0) {
        for(i in 1:n_res) {
            ## Get the columns of interest
            cells <- (1:n_traits^2)+(n_traits^2)*((n_ran) + i-1)
            ## Get the mean estimates
            template$R[[i]]$V <- matrix(apply(chain$VCV[-burnin_n, cells], 2, mean), n_traits, n_traits, byrow = FALSE)
        }
    }

    return(template)
}


