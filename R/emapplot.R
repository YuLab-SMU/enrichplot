##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "enrichResult"),
          function(x, showCategory = 30,  ...) {
              emapplot_internal(x, showCategory = showCategory, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "gseaResult"),
          function(x, showCategory = 30,  ...) {
              emapplot_internal(x, showCategory = showCategory, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 30, ...) {

              emapplot_internal(x, showCategory = showCategory, ...)
          })




##' @rdname emapplot
#' @param layout igraph layout
#' @param color Variable that used to color enriched terms, e.g. 'pvalue',
#' 'p.adjust' or 'qvalue'.
#' @param size_category relative size of the categories
#' @param min_edge The minimum similarity threshold for whether 
#' two nodes are connected, should between 0 and 1, default value is 0.2.
#' @param color_edge color of the network edge
#' @param size_edge relative size of edge width
#' @param node_label Select which labels to be displayed,
#' one of 'category', 'group', 'all' and 'none'.
#' @param pie one of 'equal' or 'Count' to set the slice ratio of the pies
#' @param group logical, if TRUE, group the category.
#' @param group_style style of ellipse, one of "ggforce" an "polygon".
#' @param label_group_style style of group label, one of "shadowtext" and "ggforce".
#' @param label_format a numeric value sets wrap length, alternatively a custom function to format axis labels.
#' @param clusterFunction function of Clustering method, such as stats::kmeans(the default),
#' cluster::clara, cluster::fanny or cluster::pam.
#' @param nWords Numeric, the number of words in the cluster tags, the default value is 4.
#' @param nCluster Numeric, the number of clusters, 
#' the default value is square root of the number of nodes.
#' @importFrom ggplot2 scale_size
#' @importFrom ggtangle geom_edge
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggrepel geom_label_repel
#' @importFrom DOSE geneInCategory
#' @author Guangchuang Yu
emapplot_internal <- function(
        x, layout=igraph::layout_with_kk, 
        showCategory = 30, color = "p.adjust",
        size_category =1,
        min_edge =.2,
        color_edge = "grey", size_edge=.5,
        node_label = "category",
        pie = "equal",
        group = FALSE,
        group_style = "ggforce",
        label_group_style = "shawdowtext",
        label_format = 30,
        clusterFunction = stats::kmeans,
        nWords = 4,
        nCluster = NULL
    ) {
    
    if (inherits(x, 'compareClusterResult')) {
        gg <- graph_from_compareClusterResult(
            x, showCategory = showCategory, color = color,
            min_edge = min_edge, size_edge = size_edge
        )
    } else {
        gg <- graph_from_enrichResult(
            x, showCategory = showCategory, color = color,
            min_edge = min_edge, size_edge = size_edge
        )        
    }    

    g <- gg$graph
    size <- vapply(gg$geneSet, length, FUN.VALUE= numeric(1))
    V(g)$size = size[V(g)$name]

    p <- ggplot(g, layout = layout) + geom_edge(color = color_edge, size = size_edge) 

    if (inherits(x, 'compareClusterResult')) {
        p <- add_node_pie(p, gg$data, pie)
    } else {
        p <- p %<+% x[, c("Description", color)] + 
            geom_point(aes(color=.data[[color]], size=.data$size)) + 
            scale_size(range=c(3, 8) * size_category) 
        p <- p + set_enrichplot_color(colors = get_enrichplot_color(2))
        p <- p + guides(size  = guide_legend(order = 1), 
                        color = guide_colorbar(order = 2, reverse = TRUE))
    }

    if (group) {
        if (inherits(x, 'compareClusterResult')) {
            p <- p + ggnewscale::new_scale_fill()
        } else {
            p <- p + ggnewscale::new_scale_color()
        }
        ggData <- groupNode(p, as.data.frame(x), nWords, clusterFunction = clusterFunction, nCluster=nCluster) 
        p$data <- ggData
        p <- add_ellipse(p, group_legend = TRUE, label_style = label_group_style, ellipse_style = group_style)
    }

    ## add node label
    if (node_label == "all" || node_label == "category") 
        p <- p + geom_text_repel(aes(label=.data$label), bg.color="white", bg.r=.1)
    ## add group label
    if (node_label == "all" || node_label == "group") {   
        label_location <- get_label_location(ggData = ggData, label_format = label_format)
        p <- p + geom_text_repel(aes(x=.data$x, y = .data$y, label=.data$label), data = label_location, bg.color="white", bg.r=.1) 
    }

    p + coord_equal() + 
        guides(size  = guide_legend(order = 1), 
               color = guide_colorbar(order = 2))
}

graph_from_enrichResult <- function(
        x, 
        showCategory = 30, color = "p.adjust",
        min_edge =.2, size_edge=.5
    ) {
    n <- update_n(x, showCategory)
    y <- as.data.frame(x)
    ## get graph.data.frame() object
    g <- get_igraph(x=x, nCategory=n, color=color, cex_line = size_edge,
                    min_edge=min_edge)
    gs <- extract_geneSets(x, n)                
    return(list(graph = g, geneSet = gs))
}

graph_from_compareClusterResult <- function(
        x,  
        showCategory = 30, color = "p.adjust",
        min_edge =.2, size_edge=.5
    ) {

    d <- tidy_compareCluster(x, showCategory)
    mergedEnrichDf <- merge_compareClusterResult(d)
    gs <- setNames(strsplit(as.character(mergedEnrichDf$geneID), "/",
                        fixed = TRUE), mergedEnrichDf$ID) 

    g <- build_emap_graph(enrichDf=mergedEnrichDf,geneSets=gs,color=color,
                        cex_line=size_edge, min_edge=min_edge, 
                        pair_sim = x@termsim, method = x@method)
    return(list(graph = g, geneSet = gs, data = d))
}



