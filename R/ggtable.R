##' plot table
##'
##'
##' @title ggtable
##' @param d data frame
##' @param p ggplot object to extract color to color rownames(d), optional
##' @return ggplot object
##' @export
##' @author guangchuang yu
ggtable <- function(d, p = NULL) {
    # has_package("ggplotify")
    ggplotify::as.ggplot(tableGrob2(d, p))
}

##' @importFrom grid gpar
##' @importFrom ggplot2 ggplot_build
tableGrob2 <- function(d, p = NULL) {
    # has_package("gridExtra")
    d <- d[order(rownames(d)),]
    tp <- gridExtra::tableGrob(d)
    if (is.null(p)) {
        return(tp)
    }

    # Fix bug: The 'group' order of lines and dots/path is different
    p_data <- ggplot_build(p)$data[[1]]
    # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
    p_data <- p_data[order(p_data[["group"]]), ]
    pcol <- unique(p_data[["colour"]])
    ## This is fine too
    ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
    j <- which(tp$layout$name == "rowhead-fg")

    for (i in seq_along(pcol)) {
        tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
    }
    return(tp)
}

