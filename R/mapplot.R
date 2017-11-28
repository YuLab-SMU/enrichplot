##' Enrichment Map
##'
##'
##' @title mapplot
##' @param x enrichment result. e.g. instance of gseaResult or enrichResult
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue, p.adjust or qvalue
##' @param layout layout of the map
##' @param ... additional parameter
##' @return ggplot object
##' @importFrom DOSE geneInCategory
##' @importFrom igraph graph.empty
##' @importFrom igraph add_vertices
##' @importFrom igraph graph.data.frame
##' @importFrom igraph delete.edges
##' @importFrom igraph V "V<-"
##' @importFrom igraph E "E<-"
##' @importFrom reshape2 melt
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' mapplot(x)
##' @author Guangchuang Yu
mapplot <- function(x, showCategory = 30, color="p.adjust", layout = "kk", ...) {
    if (!inherits(x, "gseaResult") && !inherits(x, "enrichResult"))
        stop("object not supported...")

    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    y <- y[1:n,]

    if (n == 0) {
        stop("no enriched term found...")
    } else if (n == 1) {
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- y$Description
        V(g)$color <- "red"
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    } else {
        id <- y[,1]
        geneSets <- geneSets[id]

        n <- nrow(y) #
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$Description

        for (i in 1:n) {
            for (j in i:n) {
                w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }

        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
        g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))

        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt

        colVar <- y[idx, color]
        V(g)$color <- colVar
    }


    ggraph(g, layout=layout) +
        geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey') +
        geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8))

}

