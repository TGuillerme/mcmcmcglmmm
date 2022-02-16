#' @title Plots a the roundness of a VCV matrix
#'
#' @description Displays the integral of the ordering of each major axes of an ellipse (i.e. it's roundness)
#'
#' @param VCV either a VCV list, a VCV matrix or a matrix of major axis coordinates.
#' @param dimensions which dimensions to use (if left empty, all are used).
#' @param ... Any additional arguments to be passed to plot, lines or polygon.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
## Plotting the integral
plot.roundness <- function(VCV, dimensions, ...){

    dots <- list(...)

    if(!is(VCV, "list")) {
        ## Convert into a list
        if(unique(dim(VCV)) == length(diag(VCV))) {
            VCV <- list(VCV = VCV, loc = rep(0, length(diag(VCV))), rotation = NULL)
        } else {
            VCV <- convert.VCV(VCV)
        }
    }

    ## The dimensions to check
    if(missing(dimensions)) {
        dimensions <- 1:ncol(VCV$VCV)
    }

    ## Getting the scaled major axes lengths
    majors <- sapply(dimensions, function(dim, VCV) {dist(get.one.axis(VCV, axis = dim))}, VCV = VCV)
    majors <- majors/max(majors)
    positions <- 1:length(majors)-1
    positions <- positions/max(positions)
    roundness <- sum(diff(positions)*zoo::rollmean(sort(majors), 2))

    ## Plot all that stuff
    plot_args <- list(NULL, xlim = c(0, 1), ylim = c(0, 1))
    plot_args <- get.dots(dots, plot_args, "plot", "ylab", "Scaled major axis")
    plot_args <- get.dots(dots, plot_args, "plot", "xlab", "Major axis number (scaled)")
    plot_args <- get.dots(dots, plot_args, "plot", "main", paste0("Roundness = ", round(roundness, 3)))

    ## Empty plot
    do.call(plot, plot_args)

    ## Add the polygon
    poly_args <- list(x = c(positions, rev(positions)), y = c(majors, rep(0, length(majors))))
    poly_args[c("main", "xlim", "ylim", "xlab", "ylab", "xaxt", "yaxt")] <- NULL
    poly_args <- get.dots(dots, poly_args, "polygon", "col", "grey")
    poly_args <- get.dots(dots, poly_args, "polygon", "border", NULL)
    do.call(polygon, poly_args)

    ## Add the major axes
    lines_args <- plot_args
    lines_args[c("main", "xlim", "ylim", "xlab", "ylab", "xaxt", "yaxt")] <- NULL
    lines_args <- get.dots(dots, lines_args, "lines", "col", rep("black", length(majors)))
    lines_args <- get.dots(dots, lines_args, "lines", "lwd", rep(3, length(majors)))
    lines_args <- get.dots(dots, lines_args, "lines", "lty", rep(1, length(majors)))

    for(i in 1:length(majors)) {
        one_line <- list(x = c(positions[i],positions[i]), y = c(0, majors[i]))
        one_line$lwd <- lines_args$lwd[i]
        one_line$lty <- lines_args$lty[i]
        one_line$col <- lines_args$col[i]
        do.call(lines, one_line)
    }
}