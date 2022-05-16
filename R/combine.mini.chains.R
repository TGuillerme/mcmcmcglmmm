#' @title combine.mini.chains
#'
#' @description Combine the mini chains results into a big MCMCglmm object (MCMCglmm)
#'
#' @param mini.chains a \code{"mini.chains"} objects or the path to a mini.chains object
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
combine.mini.chains <- function(mini.chains) {

    if(is(mini.chains, "character")) {
        ## Mini.chains is the path
        load(mini.chains)
        mini.chains <- output
    }

    ## Get the start end and thin values
    mcpar <- do.call(rbind, lapply(mini.chains, function(x) attr(x$Sol, "mcpar")))
    iterations_total <- sum(apply(mcpar[, c(1,2)], 1, diff))
    mcpar <- mcpar[1, ]
    mcpar[2] <- mcpar[1] + iterations_total

    ## Placeholder (with empties)
    empties <- unlist(lapply(mini.chains[[1]], is.null))
    output <- mini.chains[[1]]

    ## Get the solutions
    output$Sol <- as.mcmc(do.call(rbind, lapply(mini.chains, `[[`, "Sol")), mcpar)
    output$VCV <- as.mcmc(do.call(rbind, lapply(mini.chains, `[[`, "VCV")), mcpar)

    ## Get the deviances
    output$Deviance <- as.mcmc(do.call(c, lapply(mini.chains, `[[`, "Deviance")), mcpar)
    output$DIC <- mean(do.call(c, lapply(mini.chains, `[[`, "DIC")))

    ## Return the whole bunch
    class(output) <- class(mini.chains[[1]])
    return(output)
}


## Making x into a mcmc printable object
as.mcmc <- function(x, mcpar) {
    class(x) <- "mcmc"
    attr(x, "mcpar") <- mcpar
    return(x)
}