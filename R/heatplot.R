##' heatmap like plot for functional classification
##'
##'
##' @title heatplot
##' @param x enrichment result. e.g. instance of gseaResult or enrichResult
##' @param showCategory number of enriched terms to display
##' @param foldChange fold Change
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 theme_minimal
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_text
##' @export
##' @return ggplot object
##' @author guangchuang yu
heatplot <- function(x, showCategory=30, foldChange=NULL) {
    n <- showCategory
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    if (nrow(y) < n) {
        n <- nrow(y)
    }
    y <- y[1:n,]
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description

    readable <- x@readable
    organism <- x@organism
    if (readable & (!is.null(foldChange) ) ){
        gid <- names(foldChange)
        if (is(x, 'gseaResult')) {
            ii <- gid %in% names(x@geneList)
        } else {
            ii <- gid %in% x@gene
        }
        gid[ii] <- x@gene2Symbol[gid[ii]]
        names(foldChange) <- gid
    }

    d <- list2df(geneSets)

    if (!is.null(foldChange)) {
        d$foldChange <- fc <- foldChange[d[,2]]
        if (all(fc > 0, na.rm=TRUE)) {
            palette <- color_palette(c("blue", "red"))
        } else if (all(fc < 0, na.rm=TRUE)) {
            palette <- color_palette(c("green", "blue"))
        } else {
            palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
        }
        p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(aes_(fill=~foldChange), color="white") +
            scale_fill_gradientn(name = "fold change", colors=palette)

    } else {
        p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(color='white')
    }
    p + xlab(NULL) + ylab(NULL) + theme_minimal() +
        theme(panel.grid.major=element_blank(), axis.text.x=element_text(angle=60, hjust=1))
}
