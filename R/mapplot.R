
mapplot <- function(x, n = 50, layout = "kk", color="p.adjust", ...) {
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    if (nrow(y) < n) {
        n <- nrow(y)
    }
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
        g <- graph.data.frame(wd[,-3], directed=F)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))

        cnt <- sapply(geneSets, length)
        V(g)$size <- cnt

        colVar <- y[, color]
        V(g)$color <- colVar
    }


    ggraph(g, layout=layout) +
        geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey') +
        geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_gradientn(name = color, colors=heatmap_palette, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8))

}

heatmap_palette <- palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)


overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

