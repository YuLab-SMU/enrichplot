##' gene-concept network
##'
##'
##' @title cnetplot
##' @param x enrichment result. e.g. instance of gseaResult or enrichResult
##' @param showCategory number of enriched terms to display
##' @param foldChange fold Change
##' @param layout layout of the network
##' @param ... additional parameters
##' @return ggplot object
##' @export
##' @author guangchuang yu
cnetplot <- function(x,
                     showCategory = 5,
                     foldChange   = NULL,
                     layout = "kk",
                     ...) {

    n <- showCategory
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    if (nrow(y) < n) {
        n <- nrow(y)
    }
    y <- y[1:n,]
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description

    g <- list2graph(geneSets)

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

    size <- sapply(geneSets, length)
    V(g)$size <- min(size)/2
    V(g)$size[1:n] <- size


    if (!is.null(foldChange)) {
        fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
        V(g)$color <- NA
        V(g)$color[(n+1):length(V(g))] <- fc
        if (all(fc > 0, na.rm=TRUE)) {
            palette <- color_palette(c("blue", "red"))
        } else if (all(fc < 0, na.rm=TRUE)) {
            palette <- color_palette(c("green", "blue"))
        } else {
            palette <- color_palette(c("green", "blue", "red"))
        }
        p <- ggraph(g, layout=layout) +
            geom_edge_link(alpha=.8, colour='darkgrey') +
            geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
            scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494")
    } else {
        V(g)$color <- "#B3B3B3"
        V(g)$color[1:n] <- "#E5C494"
        p <- ggraph(g, layout=layout) +
            geom_edge_link(alpha=.8, colour='darkgrey') +
            geom_node_point(aes_(color=~I(color), size=~size))
    }

    p + scale_size(range=c(3, 8), breaks=seq(min(size), max(size), length.out=4)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void()
}


##' convert a list of gene IDs to igraph object.
##'
##'
##' @title convert gene IDs to igraph object
##' @param inputList a list of gene IDs
##' @return a igraph object.
##' @importFrom igraph graph.data.frame
##' @author Guangchuang Yu
list2graph <- function(inputList) {
    x <- data.frame()
    for (i in 1:length(inputList)) {
        x=rbind(x,
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]]))
    }
    g <- graph.data.frame(x, directed=F)
    return(g)
}


