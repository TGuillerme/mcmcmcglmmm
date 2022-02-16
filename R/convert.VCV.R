#' @title Convert VCV
#'
#' @description convert a VCV matrix into coordinates
#'
#' @param input either a VCV list, a VCV matrix or a matrix of major axis coordinates.
#' @param level which confidence interval level to use
#' 
#' @returns
#' If the input is a matrix, returns a VCV list.
#' If the input is a VCV list or a VCV matrix, returns a major axis coordinates matrix.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

## Convert a VCV to a matrix or the other way around
convert.VCV <- function(input, level = 0.95) {
    ## Detect whether the input is a matrix or a VCV

    ## Detect what to do
    do_matrix <- do_VCV <- FALSE
    if(is(input, "list") && ("VCV" %in% names(input))) {
        do_matrix <- TRUE
    } else {
        if(length(unique(dim(input))) == 1) {
            do_matrix <- TRUE
            input <- list(VCV = input, loc = rep(0, length(diag(VCV))))
        } else {
            do_VCV <- TRUE
        }
    }

    ## Make a matrix
    if(do_matrix) {
        ## Check if location
        if(is.null(input$loc)) {
            if(length(rename <- which(names(input) %in% "Sol")) > 0) {
                ## Change Sol to loc
                names(input)[rename] <- "loc"
            } else {
                input$loc <- rep(0, diag(input$VCV))
            }
        }
    
        ## Get the major axis
        matrix <- get.one.axis(input, axis = c(1:length(diag(input$VCV))), level = level)
        ## Add the matrix row names
        rownames(matrix) <- 1:nrow(matrix)
        ## Name them correctly
        rownames(matrix)[1:(nrow(matrix)/2)] <- paste0("axis", 1:(nrow(matrix)/2), ".1")
        rownames(matrix)[((nrow(matrix)/2)+1):nrow(matrix)] <- paste0("axis", 1:(nrow(matrix)/2), ".2")
        colnames(matrix) <- paste0("D", 1:ncol(matrix))
        ## Rotate the matrix
        if(!is.null(input$rotation)) {
            matrix <- matrix %*% input$rotation
        }
        ## Return the matrix
        return(matrix)
        # plot(ellipse::ellipse(VCV), type = "l", xlim = c(-3, 3), ylim =c(-3, 3))
        # lines(matrix[c("axis1.1","axis1.2"), ])
        # lines(matrix[c("axis2.1","axis2.2"), ])
    }

    if(do_VCV) {
        ## Reverse the process (might be bugged!)
        return(var(input/sqrt(qchisq(level, df = 1))))
    }
    stop("Input must be a VCV matrix, a VCV list or a matrix of coordinates.")
}
