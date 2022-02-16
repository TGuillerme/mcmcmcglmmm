#' @title Making flat priors
#'
#' @description Making a list of flat priors to be fed to a mini.chain
#'
#' @param residuals the number of residual terms (R structure in the list)
#' @param randoms the number of random terms (G structure in the list)
#' @param ntraits the number of traits
#' @param nu the value for nu (applied throughout)
#' 
#' @examples
#' ## Default flat priors for 2 traits
#' flat.prior(ntraits = 2)
#' 
#' ## Flat priors for a model with 1 residual and 4 random terms
#' ## for 3 traits and a nu of 0.002
#' flat.prior(residuals = 1, randoms = 4, ntraits = 3, nu = 0.002)
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

flat.prior <- function(residuals = 1, randoms = 1, ntraits, nu = 0) {

    ## Making the prior template
    flat_prior_template <- list(V = diag(ntraits), nu = nu)

    ## Making the residuals list
    if(residuals != 0) {
        R <- list()
        for(i in 1:residuals) {
            R <- c(R, list(flat_prior_template))
            names(R)[i] <- paste0("R", i)
        }
    } else {
        R <- NULL
    }
    ## Making the randoms list
    if(randoms != 0) {
        G <- list()
        for(i in 1:randoms) {
            G <- c(G, list(flat_prior_template))
            names(G)[i] <- paste0("G", i)
        }
    } else {
        G <- NULL
    }

    ## Return the flat priors
    return(list(R = R, G = G))
}