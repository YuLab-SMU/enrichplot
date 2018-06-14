##' plot table
##'
##'
##' @title ggtable
##' @param d data frame
##' @param p ggplot object to extract color to color rownames(d), optional
##' @return ggplot object
##' @importFrom ggplotify as.ggplot
##' @export
##' @author guangchuang yu
ggtable <- function(d, p = NULL) {
    as.ggplot(tableGrob2(d, p))
}

##' @importFrom grid gpar
##' @importFrom gridExtra tableGrob
##' @importFrom ggplot2 ggplot_build
tableGrob2 <- function(d, p = NULL) {
    d <- d[order(rownames(d)),]
    tp <- tableGrob(d)
    if (is.null(p)) {
        return(tp)
    }
    pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
    j <- which(tp$layout$name == "rowhead-fg")

    for (i in seq_along(pcol)) {
        tp$grobs[j][[i+1]][["gp"]] = gpar(col = pcol[i])
    }
    return(tp)
}

