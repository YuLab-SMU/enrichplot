##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "enrichResult"),
          function(x, showCategory = 5,
                   foldChange = NULL, layout = "kk", ...) {
              cnetplot.enrichResult(x, showCategory = showCategory,
                                    foldChange = foldChange, layout = layout, ...)
          })

##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "gseaResult"),
          function(x, showCategory = 5,
                   foldChange = NULL, layout = "kk", ...) {
              cnetplot.enrichResult(x, showCategory = showCategory,
                                    foldChange = foldChange, layout = layout, ...)
          })

##' @rdname cnetplot
##' @param colorEdge whether coloring edge by enriched terms
##' @param circular whether using circular layout
##' @param node_label whether display node label
##' @importFrom ggraph geom_edge_arc
##' @author Guangchuang Yu
cnetplot.enrichResult <- function(x,
                     showCategory = 5,
                     foldChange   = NULL,
                     layout = "kk",
                     colorEdge = FALSE,
                     circular = FALSE,
                     node_label = TRUE,
                     ...) {

    if (circular) {
        layout <- "linear"
        geom_edge <- geom_edge_arc
    } else {
        geom_edge <- geom_edge_link
    }

    geneSets <- extract_geneSets(x, showCategory)

    g <- list2graph(geneSets)

    foldChange <- fc_readable(x, foldChange)

    size <- sapply(geneSets, length)
    V(g)$size <- min(size)/2

    n <- length(geneSets)
    V(g)$size[1:n] <- size

    if (colorEdge) {
        E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
        edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)
    } else {
        edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
    }

    if (!is.null(foldChange)) {
        fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
        V(g)$color <- NA
        V(g)$color[(n+1):length(V(g))] <- fc
        palette <- fc_palette(fc)
        p <- ggraph(g, layout=layout, circular = circular) +
            edge_layer +
            geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
            scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494")
    } else {
        V(g)$color <- "#B3B3B3"
        V(g)$color[1:n] <- "#E5C494"
        p <- ggraph(g, layout=layout, circular=circular) +
            edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size))
    }

    p <- p + scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
        theme_void()

    if (node_label)
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE)

    return(p)
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
    x <- list2df(inputList)
    g <- graph.data.frame(x, directed=FALSE)
    return(g)
}


list2df <- function(inputList) {
    ldf <- lapply(1:length(inputList), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}


