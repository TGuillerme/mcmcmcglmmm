context("make.ellipse")

## Test
test_that("make.ellipse works", {

    ## All the lambda/shape business

    ## Visual examples
    ## Run the example with 100D
    input_values <- seq(from = 0, to = 1, length.out = 11)
    shape_list <- get.lambda(input_values)
    plot(NULL, xlim = c(0,1), ylim = c(0,1), xlab = "Dimensions (proportional)", ylab = "Axis length (proportional)", main = "Axis lengths with in 100D")
    cols <- rainbow(length(shape_list))
    for(i in 1:length(shape_list)) {
        lines(generate.shape(dimensions = 100, lambda = shape_list[i]), col = cols[i])
    }
    legend(x = 0.5, y = 0.75, legend = paste(paste0("shape = ", round(input_values, 1)), paste0("lambda = ", round(shape_list, 2)), sep = "; "), lty = 1, col = cols, cex = 0.7)

    ## Run the example with 3D and minimum = 0.1
    input_values <- seq(from = 0, to = 1, length.out = 11)
    shape_list <- get.lambda(input_values)
    plot(NULL, xlim = c(0,1), ylim = c(0,1), xlab = "Dimensions (proportional)", ylab = "Axis length (proportional)", main = "Axis lengths with in 3D\nwith min.thick = 0.01")
    cols <- rainbow(length(shape_list))
    for(i in 1:length(shape_list)) {
        lines(generate.shape(dimensions = 3, lambda = shape_list[i], min.thick = 0.01), col = cols[i])
    }
    legend(x = 0.5, y = 0.75, legend = paste(paste0("shape = ", round(input_values, 1)), paste0("lambda = ", round(shape_list, 2)), sep = "; "), lty = 1, col = cols, cex = 0.7)


    ## Testing make.VCV
    ## Default is a circle
    test <- make.VCV()
    expect_equal(names(test), c("VCV", "loc", "rotation"))
    expect_null(plot.VCV(test, VCV.vectors = TRUE, eigen.vectors = TRUE, eigen.values = TRUE))

    ## Example with options
    expect_null(plot.VCV(make.VCV(shape = 0.7, covariance = 0.5), VCV.vectors = TRUE, eigen.vectors = TRUE, eigen.values = TRUE, col.eigen.vectors = "red", lty.VCV.vectors = 3, lwd = 2))
})
