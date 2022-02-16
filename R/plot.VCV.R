#' @title Plots an ellipse
#'
#' @description Plotting an ellipse from the output of make.ellipse
#'
#' @param VCV either a VCV list, a VCV matrix or a matrix of major axis coordinates.
#' @param ellipse logical, whether to plot the confidence interval ellipse or not (default is TRUE).
#' @param major.axes logical, whether to plot the major axes or not (default is TRUE).
#' @param level the confidence interval level to estimate the ellipse and the major axes (default is 0.95)
#' @param dimensions which dimensions to plot (up to 3!). Default is c(1,2).
#' @param VCV.vectors display the VCV matrix vectors (default is FALSE)
#' @param eigen.vectors display the eigen vectors (default is TRUE)
#' @param eigen.values if eigen.vectors = TRUE, whether to scale the eigen vectors by the eigen values.
#' @param add logical, whether to add it to an existing plot (TRUE) or not (FALSE; default)
#' @param ... Any additional arguments to be passed to plot, lines or arrows.
#' 
#' @details
#' ... can take specific options for specific elements to plot if you precede the plot argument by the element name. For example \code{col.eigen.vectors = "blue"} will make only the eigen.vectors blue.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

plot.VCV <- function(VCV, ellipse = TRUE, major.axes = TRUE, level = 0.95, dimensions = c(1,2), VCV.vectors = FALSE, eigen.vectors = FALSE, eigen.values = FALSE, add = FALSE, ...) {

    dots <- list(...)

    if(!is(VCV, "list")) {
        ## Convert into a list
        if(unique(dim(VCV)) == length(diag(VCV))) {
            VCV <- list(VCV = VCV, loc = rep(0, length(diag(VCV))), rotation = NULL)
        } else {
            VCV <- convert.VCV(VCV)
        }
    }

    ## Get the dimensions to plot
    is_3D <- length(dimensions) == 3

    ## Plot the ellipse
    if(is_3D) {
        ## Make the ellipse
        plot_ellipse <- rgl::ellipse3d(VCV$VCV, centre = VCV$loc, level = level)
    } else {
        ## Make the ellipse
        plot_ellipse <- reloc(ellipse::ellipse(VCV$VCV, level = level), VCV)[, dimensions]
    }

    ## Get the plot arguments
    if(ellipse) {
        plot_args <- list(x = plot_ellipse)
    } else {
        plot_args <- list(x = NULL)
    }
    if(is_3D) {
        ## Set the plotting function
        plot.fun <- rgl::plot3d
        lines.fun <- rgl::plot3d
        ## Set the parameters for the 3D ellipse
        plot_args <- get.dots(dots, plot_args, "ellipse", "col", "grey")
        plot_args <- get.dots(dots, plot_args, "ellipse", "alpha", 0.5)
    } else {
        ## Set the plotting fun
        plot.fun <- plot
        lines.fun <- lines
        ## Parameters for the 2D ellipse
        plot_args <- get.dots(dots, plot_args, "ellipse", "col", "black")
        plot_args <- get.dots(dots, plot_args, "ellipse", "lwd", 1)
        plot_args <- get.dots(dots, plot_args, "ellipse", "lty", 1)
    }


    if(!add) {
        ## Parameters for initialising the plot
        if(!is_3D) {
            plot_args <- get.dots(dots, plot_args, "ellipse", "type", "l")
        }
        plot_args <- get.dots(dots, plot_args, "ellipse", "main")

        ## Plot limits
        if(is_3D) {
            lims <- range(pretty(range(plot_ellipse$vb)))
            ## Z axis
            plot_args <- get.dots(dots, plot_args, "ellipse", "zlim", lims)
            plot_args <- get.dots(dots, plot_args, "ellipse", "zlab", paste0("Dim.", dimensions[3]))
        } else {
            lims <- range(pretty(plot_ellipse[, dimensions]))
        }
        plot_args <- get.dots(dots, plot_args, "ellipse", "xlim", lims)
        plot_args <- get.dots(dots, plot_args, "ellipse", "ylim", lims)
        plot_args <- get.dots(dots, plot_args, "ellipse", "xlab", paste0("Dim.", dimensions[1]))
        plot_args <- get.dots(dots, plot_args, "ellipse", "ylab", paste0("Dim.", dimensions[2]))

        ## Plot the ellipse
        do.call(plot.fun, plot_args)
    } else {
        if(is_3D) {
            plot_args$add <- TRUE
        }
        do.call(lines.fun, plot_args)
    }

    ## Plot the major axis
    if(major.axes) {
        ## Calculate the major axes
        major_axes <- reloc(get.one.axis(VCV$VCV, level = level, axis = dimensions), VCV)
        ## Get the plot arguments
        axes_args <- list()
        axes_args <- get.dots(dots, axes_args, "major.axes", "col", rep("black", length(dimensions)))
        axes_args <- get.dots(dots, axes_args, "major.axes", "lty", 2)
        axes_args <- get.dots(dots, axes_args, "major.axes", "lwd", 1)
        if(is_3D) {
            lines.fun <- rgl::lines3d
        } else {
            lines.fun <- lines
        }
        ## Plot each major axis
        for(i in dimensions) {
            axis_args <- axes_args
            axis_args$x <- major_axes[c(i, i+length(dimensions)), dimensions]
            axis_args$col <- axes_args$col[i]
            do.call(lines.fun, axis_args)
        }
    }


    ## Plot the VCV vectors
    if(VCV.vectors) {
        ## Get the VCV vectors
        vcv_vec <- reloc(vectors.VCV(VCV$VCV), VCV)

        ## Get the arguments
        vcv_args <- list()
        vcv_args <- get.dots(dots, vcv_args, "VCV.vectors", "col", "orange")
        vcv_args <- get.dots(dots, vcv_args, "VCV.vectors", "lwd", 1)
        vcv_args <- get.dots(dots, vcv_args, "VCV.vectors", "lty", 1)
        if(!is_3D) {
            vcv_args <- get.dots(dots, vcv_args, "VCV.vectors", "length", 0.1)
            vcv_args <- get.dots(dots, vcv_args, "VCV.vectors", "angle", 30)
            vcv_args <- get.dots(dots, vcv_args, "VCV.vectors", "code", 2)
        }
    
        ## Plot each VCV vector
        for(i in dimensions) {
            one_vcv <- vcv_args
            if(is_3D) {
                one_vcv$x <- rbind(rep(0, length(dimensions)), vcv_vec[i, ])
                do.call(segments3d, one_vcv)
            } else {
                one_vcv$x0 <- VCV$loc[dimensions[1]]
                one_vcv$y0 <- VCV$loc[dimensions[2]]
                one_vcv$x1 <- vcv_vec[i, 1]
                one_vcv$y1 <- vcv_vec[i, 2]
                do.call(arrows, one_vcv)
            }
        }
    }

    ## Plot the eigen vectors
    if(eigen.vectors) {
        ## Get the eigen vectors
        eig <- eigen(VCV$VCV)

        ## Rescale them by the values
        if(eigen.values) {
            scalar <- eig$values
        } else {
            scalar <- rep(1, length(eig$values))
        }

        ## Relocate them
        eig <- reloc(eig$vectors, VCV)
        ## Get the arguments
        eig_args <- list()
        eig_args <- get.dots(dots, eig_args, "eigen.vectors", "col", "blue")
        eig_args <- get.dots(dots, eig_args, "eigen.vectors", "lwd", 1)
        eig_args <- get.dots(dots, eig_args, "eigen.vectors", "lty", 1)
        if(!is_3D) {
            eig_args <- get.dots(dots, eig_args, "eigen.vectors", "length", 0.1)
            eig_args <- get.dots(dots, eig_args, "eigen.vectors", "angle", 30)
            eig_args <- get.dots(dots, eig_args, "eigen.vectors", "code", 2)
        }
        ## Plot each eigen vector
        for(i in dimensions) {
            one_eig <- eig_args
            if(is_3D) {
                one_eig$x <- rbind(rep(0, length(dimensions)), eig[, i])
                do.call(segments3d, one_eig)
            } else {
                one_eig$x0 <- loc[dimensions[1]]
                one_eig$y0 <- loc[dimensions[2]]
                one_eig$x1 <- eig[1, i] * scalar[i]
                one_eig$y1 <- eig[2, i] * scalar[i]
                do.call(arrows, one_eig)
            }
        }
    }
    return(invisible())
}

## Get options
get.dots <- function(dots, args, suffix, name, default) {

    if(!missing(suffix) && !is.null(names(dots))) {
        name_plus <- paste0(name, ".", suffix)
        ## Override the name without suffix in dots
        if(name_plus %in% names(dots)) {
            dots[[name]] <- dots[[name_plus]]
        }
    } 

    if(!missing(default)) {
        if(is.null(dots[[name]])) {
            args[[name]] <- default
        } else {
            args[[name]] <- dots[[name]]
        }


    } else {
        if(!is.null(dots[[name]])) {
            args[[name]] <- dots[[name]]
        }
    }
    return(args)
}
## Relocate points
reloc <- function(matrix, VCV) {
    if(!is.null(VCV$loc)) {
        ## Relocate the ellipse
        for(i in 1:ncol(matrix)) {
            matrix[,i] <- matrix[,i] + VCV$loc[i]
        }
    }
    return(matrix)
}