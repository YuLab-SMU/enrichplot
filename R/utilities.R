
##' create color palette for continuous data
##'
##'
##' @title color_palette
##' @param colors colors of length >=2
##' @return color vector
##' @importFrom grDevices colorRampPalette
##' @export
##' @author guangchuang yu
color_palette <- function(colors) colorRampPalette(colors)(n = 299)

heatmap_palette <- color_palette(c("red", "yellow", "green"))

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}



