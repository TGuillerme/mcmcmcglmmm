#' @title Get lambda
#'
#' @description Get the lambda value corresponding to the input shape
#'
#' @param shape the shape value
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

## Internals
## This function approximates the lambda value corresponding the the correct shape
## Shape = 0: circle          -> lambda = Inf
## Shape = (0;0.5): pancake   -> lambda = (Inf; 1)
## Shape = 0.5: clean ellipse -> lambda = 1
## Shape = (0.5;1): cigar     -> lambda = (1; 0)
get.lambda <- function(shape) {
    ## Reverse x
    shape <- 1-shape
    ## Transform the extremes
    shape <- ifelse(shape == 0, -Inf, shape)
    shape <- ifelse(shape == 1, Inf, shape)
    ## Return the centred exp val of shape (scaled and centred so that 0.5 = ellipse)
    return(exp((shape*4-2)*2))
}