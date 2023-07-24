#' @title diagnose.mini.chains
#'
#' @description Diagnose a set of mini chains results
#'
#' @param mini.chains a mini chain MCMCglmm result
#' @param what if mini.chains is a single chain, which part from \code{\link[MCMCglmm]{summary.MCMCglmm}} to diagnose (default is \code{"eff.samp"} for the effective sample size).
#' @param plot if mini.chains is a single chain, logical, whether to plot the results (\code{TRUE}; default) or not (\code{FALSE}).
#' @param ... if mini.chains is a single chain: if \code{plot = TRUE}, any optional arguments to be past to \code{\link[graphics]{hist}}. Else, any optional arguments to be passed to \code{\link[coda]{gelman.diag}}.
#' 
#' @details
#' If the input is a list of chains straight from \code{run.mini.chains}, a Gelman and Rubin convergence diagnosis is run (\code{\link[coda]{gelman.diag}}). Else if the input is a single chain combined from \code{combine.mini.chains} only the effective sample size for each parameter (or \code{what}) is measured.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

diagnose.mini.chains <- function(mini.chains, what = "eff.samp", plot = TRUE, ...) {

    ## Toggle the difference diagnosis options
    mini.chains_class <- class(mini.chains)
    convergence <- FALSE
    switch(mini.chains_class,
        "mini.chains" = {convergence <- length(mini.chains) > 1},
        ## Is a single MCMCglmm object
        "MCMCglmm"    = {convergence <- FALSE},
        ## Is a single mcmc object
        "mcmc"        = {mini.chains <- list(mini.chains); convergence <- FALSE})
    if(mini.chains_class == "mini.chains" && !convergence) {
        mini.chains <- mini.chains[[1]]
    }

    ## Do the convergence diagnosis
    if(convergence) {
        args <- list(x = lapply(mini.chains, `[[`, 1), ...)
        return(do.call(coda::gelman.diag, args))
    }

    ## Summarising the chain
    summarised <- summary.MCMCglmm(mini.chains)

    ## How many bits are summarised?
    # elements <- names(summarised)
    elements <- c("solutions", "Gcovariances", "Rcovariances") # TODO: make that automatic
    diagnoses <- lapply(as.list(elements), function(one_element, summarised, what) return(summarised[[one_element]][, what]), summarised, what)
    names(diagnoses) <- elements
    ## Remove the empties
    if(any(nulls <- unlist(lapply(diagnoses, is.null)))) {
        diagnoses <- diagnoses[-nulls]
        elements <- elements[-nulls]
    }

    ## Plot the diagnosis
    if(plot) {
        ## Setting the plot sizes
        op <- par(mfrow = c(length(diagnoses),1))

        ## Do one plot
        one.plot <- function(ESS, main, what, ...) {
            graphics::hist(ESS, main = main, ...)
            if(what == "eff.samp") {
                graphics::abline(v = 200)
            }
        }

        ## Do all the plots
        silent <- mapply(one.plot, diagnoses, elements, MoreArgs = list(what = what, ...), SIMPLIFY = FALSE)
        par(op)
    }

    return(diagnoses)
}