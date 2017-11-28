
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

sig_palette <- color_palette(c("red", "yellow", "blue"))

heatmap_palette <- color_palette(c("red", "yellow", "green"))

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

fc_readable <- function(x, foldChange = NULL) {
    if (is.null(foldChange))
        return(NULL)

    if(x@readable) {
        gid <- names(foldChange)
        if (is(x, 'gseaResult')) {
            ii <- gid %in% names(x@geneList)
        } else {
            ii <- gid %in% x@gene
        }
        gid[ii] <- x@gene2Symbol[gid[ii]]
        names(foldChange) <- gid
    }
    return(foldChange)
}

fc_palette <- function(fc) {
    if (all(fc > 0, na.rm=TRUE)) {
        palette <- color_palette(c("blue", "red"))
    } else if (all(fc < 0, na.rm=TRUE)) {
        palette <- color_palette(c("green", "blue"))
    } else {
        palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
    }
    return(palette)
}

update_n <- function(x, showCategory) {
    n <- showCategory
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    if (nrow(x) < n) {
        n <- nrow(x)
    }
    return(n)
}

extract_geneSets <- function(x, n) {
    n <- update_n(x, n)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    y <- y[1:n,]
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description
    return(geneSets)
}

