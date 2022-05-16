#' @title make.mini.chains
#'
#' @description Making a list of mini.chains object (MCMCglmm primers)
#'
#' @param data a full ordinated matrix with a column called \code{"animal"} and potentially one called \code{"clade"} (if needed).
#' @param dimensions a numeric or character vector of the dimensions to include.
#' @param tree the \code{"phylo"} or \code{"multiPhylo"} object to run the models.
#' @param trait.family the family of the traits (default is \code{"gaussian"}).
#' @param residuals optional, the type of residuals (see details).
#' @param randoms optional, the type of randoms (see details).
#' @param parameters the named list of parameters for the MCMCglmm (if missing, it is set to \code{list(nitt=1000, burnin = 100, thin = 100)} for a fast test).
#' @param priors A list of priors or a value of the overal nu parameter for generating flat priors for the \code{\link[MCMCglmm]{MCMCglmm}}.
#' @param verbose whether to make the \code{\link[MCMCglmm]{MCMCglmm}} verbose (\code{TRUE}; default) or not (\code{FALSE}).
#' @param ... any other parameters to be passed to \code{\link[MCMCglmm]{MCMCglmm}}.
#' 
#' @details
#' By default, no model is used for the random terms (\code{NULL}). The types of model for the residuals and random terms can be a combination of any of the following:
#' \itemize{
#'  \item \code{"global"} for just an overall term effect. For the random terms that would be just a phylogenetic effect (\code{~us(trait):animal}) or for the residuals just an overall category (\code{~us(trait):unit}).
#'  \item the name of the column(s) containing the factors for just as many levels of effects as there are non-empty levels in the data (\code{""} or \code{<NA>} are considered as empty levels). For example for the column \code{"clade"} the added random effects will be something like \code{us(at.level(clade,1):trait):animal} (i.e. a special phylogenetic random effect for clade 1) or \code{us(at.level(clade,1):trait):unit} for the residuals (i.e. a special residual effect for clade 1).
#' }
#' 
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
make.mini.chains <- function(data, dimensions, tree, trait.family = "gaussian", residuals = "global", randoms = NULL, parameters, priors = 0.02, verbose = TRUE, ...) {

    ## Setting the fixed effect model
    ## Getting the terms
    rand_terms <- handle.terms(randoms)
    resi_terms <- handle.terms(residuals)
    if(is.null(resi_terms$main) && is.null(resi_terms$add)) {
        stop("needs a residual term")
    }

    ## Is there any clade effect in the model?
    if(is.null(resi_terms$add) && is.null(rand_terms$add)) {
        other_terms <- FALSE
        ## Setting the fixed model
        fixed_model <- ~trait-1
    } else {
        other_terms <- TRUE
        if(!is.null(resi_terms$add) && !is.null(rand_terms$add)) {
            adding_terms <- unique(randoms[rand_terms$add], residuals[resi_terms$add])
        } else {
            if(!is.null(rand_terms$add)) {
                adding_terms <- unique(randoms[rand_terms$add])
            } else {
                adding_terms <- unique(residuals[resi_terms$add])
            }
        }

        ## Setting the fixed model
        fixed_model <- as.formula(paste("~", paste(paste0("trait:", adding_terms), collapse = " + "), "- 1"))
    }

    ## Find the number of clades
    if(other_terms) {
        if(any(wrongs <- !(adding_terms %in% colnames(data)))) {
            stop(paste0("The following column", ifelse(sum(wrongs) > 1, "s ", " "), "cannot be found in data:", paste(adding_terms[wrongs], collapse = ", "), "."))
        } 
        ## Check if the other terms are factors
        n_groups <- list()
        for(i in 1:length(adding_terms)) {
            if(!is(data[, adding_terms[i]], "factor")) {
                data[, adding_terms[i]] <- as.factor(data[, adding_terms[i]])     
            }
            ## Getting the group numbers
            n_groups[[i]] <- seq_along(levels(data[, adding_terms[i]]))
            ## Removing empties
            empties <- which(levels(data[, adding_terms[i]]) == "" | is.na(levels(data[, adding_terms[i]])))
            if(length(empties) > 0) {
                n_groups[[i]] <- n_groups[[i]][-empties]
            }
            ## Making the level names
            n_groups[[i]] <- paste(adding_terms[i], n_groups[[i]], sep = ",")
        }
        n_groups <- unname(unlist(n_groups))
    }

    ## Check the dimensionality of the data
    if(!all(dim_check <- apply(data[, dimensions, drop = FALSE], 2, is, "numeric"))) {
        stop(paste0("Invalid dimensions (not numeric?): ", paste(names(which(!dim_check)), collapse = ", "), "."))
    }

    ## Change the dimensions to numeric
    if((class_dim <- class(dimensions)) != "numeric") {
        select_dim <- switch(class_dim,
                             "character" = match(dimensions, colnames(data)),
                             "integer"   = as.numeric(dimensions))
    } else {
        select_dim <- dimensions
    }
    ## Reduce the dataset size
    non_data <- unlist(lapply(as.list(1:ncol(data)), function(col, data) return(class(data[,col])), data = data))
    data_reduce <- data[, c(select_dim, which(non_data != "numeric")), drop = FALSE]
    
    ## Match the data and the trees
    cleaning <- clean.data(data_reduce, tree)
    data_reduce <- cleaning$data
    tree <- cleaning$tree
    ## Tell what happens
    if(any(!is.na(cleaning$dropped_tips))) {
        warning(paste0("Dropped ", length(cleaning$dropped_tips), ifelse(length(cleaning$dropped_tips) == 1, " tip", " tips"), " from the ", ifelse(is(tree, "multiPhylo"), "trees", "tree"), " that ", ifelse(length(cleaning$dropped_tips) == 1, "was", "were"), " not present in the data."))
    }
    if(any(!is.na(cleaning$dropped_rows))) {
        warning(paste0("Dropped ", length(cleaning$dropped_rows), ifelse(length(cleaning$dropped_rows) == 1, " row", " rows"), " from the dataset that ", ifelse(length(cleaning$dropped_rows) == 1, "was", "were"), " not present in the ", ifelse(is(tree, "multiPhylo"), "trees", "tree") , "."))
    }

    ## Fixed effect
    ## Setting the fixed formula (initialising)
    fixed <- fixed_model
    ## Moving the model to the second component of the formula (third)
    fixed[[length(fixed_model) + 1]] <- fixed[[length(fixed_model)]]
    ## Adding the trait model as the first component (first)
    fixed[[length(fixed_model)]] <- as.formula(paste0("~cbind(", paste(colnames(data)[select_dim], collapse = ", "), ")"))[[2]]

    ## Random effect
    random <- NULL
    n_randoms <- 0
    if(!is.null(rand_terms$main)) {
        if(!is.null(rand_terms$add)) {
            ## Group terms + global
            random <- group.terms(n_groups, type = "animal", add.global = TRUE)
            n_randoms <- length(n_groups) + 1            
        } else {
            ## Global terms only
            random <- ~ us(trait):animal        
            n_randoms <- 1
        }
    } else {
        if(!is.null(rand_terms$add)) {
            ## Group terms only
            random <- group.terms(n_groups, type = "animal", add.global = FALSE)
            n_randoms <- length(n_groups)         
        }   
    }

    ## Residuals effect
    rcov <- NULL
    n_residuals <- 0
    if(!is.null(resi_terms$main)) {
        if(!is.null(resi_terms$add)) {
            ## Group terms + global
            rcov <- group.terms(n_groups, type = "units", add.global = TRUE)
            n_residuals <- length(n_groups) + 1            
        } else {
            ## Global terms only
            rcov <- ~ us(trait):units       
            n_residuals <- 1
        }
    } else {
        if(!is.null(resi_terms$add)) {
            ## Group terms only
            rcov <- group.terms(n_groups, type = "units", add.global = FALSE)
            n_residuals <- length(n_groups)         
        }   
    }

    ## Priors
    if(is(priors, "numeric")) {
        priors <- flat.prior(ntraits = length(select_dim), residuals = n_residuals, randoms = n_randoms, nu = priors)
    }

    ## Parameters
    if(missing(parameters)) {
        parameters <- list()
    }
    if(is.null(parameters$nitt)) {
        parameters$nitt <- 100
    }
    if(is.null(parameters$burnin)) {
        parameters$burnin <- 10
    }
    if(is.null(parameters$thin)) {
        parameters$thin <- 10
    }    

    ## Distributions
    if(length(trait.family) == 1) {
        family <- rep(trait.family, length(select_dim))
    } else {
        if(length(trait.family) == length(select_dim)) {
            family <- trait.family
        } else {
            stop("Incorrect family (must be either a single character string or of the same number as dimensions).")
        }
    }

    ## Set the tree
    if(is(tree, "phylo")) {
        tree <- list(tree)
        class(tree) <- "multiPhylo"
    }
    if(!(is(tree, "multiPhylo"))) {
        stop("tree must be a phylo or multiPhylo object.")
    }

    ## Setting the tree(s)
    params_list <- list(fixed = fixed, random = random, rcov = rcov, family = family, priors = priors, verbose = verbose, parameters = parameters, ...)
    output <- list(data   = data_reduce,
                   tree   = tree,
                   params = params_list,
                   run    = function(one_data, one_tree, params) MCMCglmm(
                            fixed    = params$fixed,
                            random   = params$random,
                            rcov     = params$rcov,
                            family   = params$family,
                            pedigree = one_tree,
                            data     = one_data,
                            prior    = params$priors,
                            verbose  = params$verbose,
                            burnin   = params$parameters$burnin,
                            nitt     = params$parameters$nitt,
                            thin     = params$parameters$thin,
                            saveX    = FALSE,
                            saveZ    = FALSE,
                            saveXL   = FALSE,
                            pl       = FALSE,
                            ...)
                 )
    class(output) <- c("mini.chains")
    return(output)
}

## Making the group terms 
group.terms <- function(n_groups, type, add.global = FALSE) {
    
    ## Term for the first clade
    form <- paste0("~us(at.level(", n_groups[1], "):trait):", type)

    ## Term for the subsequent clades
    for(i in 2:length(n_groups)) {
        form <- paste0(form, paste0(" +us(at.level(", n_groups[i], "):trait):", type))
    }
    ## Add the global term?
    if(add.global) {
        add_global <- 1
        form <- paste0(form, paste0(" +us(trait):", type))
    } else {
        add_global <- 0
    }
    return(as.formula(form))
}

## Handling the terms
handle.terms <- function(term) {
    ## No terms
    if(is.null(term)) {
        return(list(main = NULL, add = NULL))
    }

    ## Get the main term (global)
    if(any("global" %in% term)) {
        main <- which(term == "global")
    } else {
        main <- NULL
    }
    ## Get the additional terms
    if(any(term != "global")) {
        add <- which(term != "global")
    } else {
        add <- NULL
    }
    ## Return the terms IDs
    return(list(main = main, add = add))
}

