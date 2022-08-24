##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "enrichResult"),
          function(x, showCategory = 30,  ...) {
              emapplot.enrichResult(x, showCategory = showCategory, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "gseaResult"),
          function(x, showCategory = 30,  ...) {
              emapplot.enrichResult(x, showCategory = showCategory, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 30, ...) {

              emapplot.compareClusterResult(x, showCategory = showCategory, ...)
          })




##' @rdname emapplot
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
##' @importFrom DOSE geneInCategory
##' @importFrom ggplot2 scale_color_discrete
##' @importFrom ggplot2 scale_size_continuous
##' @importFrom ggplot2 scale_fill_discrete
##' @importFrom ggplot2 geom_text
##' @importFrom shadowtext geom_shadowtext
##' @param layout Layout of the map, e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
##' 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
##' @param coords a data.frame with two columns: 'x' for X-axis coordinate and 'y' for Y-axis coordinate.
##' @param color Variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @param cex_line Scale of line width
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##' @param cex_label_category Scale of category node label size.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @param shadowtext a logical value, whether to use shadow font.
##' @param label_style style of group label, one of "shadowtext" and "ggforce".
##' @param repel whether to correct the position of the label. Defaults to FALSE.
##' @param node_label Select which labels to be displayed,
##' one of 'category', 'group', 'all' and 'none'.
##' @param with_edge Logical, if TRUE, draw the edges of the network diagram.
##' @param group_category a logical, if TRUE, group the category.
##' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
##' The default is FALSE.
##' @param cex_label_group Numeric, scale of group labels size, the default value is 1.
##' @param nWords Numeric, the number of words in the cluster tags, the default value is 4.
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param clusterFunction function of Clustering method, such as stats::kmeans(the default),
##' cluster::clara, cluster::fanny or cluster::pam.
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' @param ... additional parameters
##' 
##' additional parameters can refer the following parameters.
##'     \itemize{
##'        \item \code{force} Force of repulsion between overlapping text labels. Defaults to 1. 
##'        \item \code{nudge_x, nudge_y} Horizontal and vertical adjustments to nudge 
##'         the starting position of each text label. 
##'        \item \code{direction} "both", "x", or "y" – direction in which to adjust position of labels.
##'        \item \code{ellipse_style} style of ellipse, one of "ggforce" an "polygon".
##'        \item \code{ellipse_pro} numeric indicating confidence value for the ellipses, 
##'         it can be used only when ellipse_style = "polygon".
##'        \item \code{alpha} the transparency of ellipse fill.
##'        \item \code{type} The type of ellipse. The default "t" assumes a multivariate t-distribution, 
##'         and "norm" assumes a multivariate normal distribution. "euclid" draws a circle with the 
##'         radius equal to level, representing the euclidean distance from the center. 
##'     }
##' 
##' @author Guangchuang Yu
emapplot.enrichResult <- function(x, showCategory = 30, 
                                  layout = NULL, coords = NULL,
                                  color = "p.adjust", min_edge = 0.2,
                                  cex_label_category  = 1, cex_category = 1,
                                  cex_line = 1, shadowtext = TRUE,
                                  label_style = "shadowtext", repel = FALSE,
                                  node_label  = "category",
                                  with_edge = TRUE, group_category = FALSE,  
                                  group_legend = FALSE,                             
                                  cex_label_group = 1, nWords = 4, 
                                  label_format = 30,
                                  clusterFunction = stats::kmeans,
                                  nCluster = NULL,  ...) {
    has_pairsim(x)
    label_size_category <- 5
    label_group <- 3
    n <- update_n(x, showCategory)
    y <- as.data.frame(x)
    ## get graph.data.frame() object
    g <- get_igraph(x=x, nCategory=n, color=color, cex_line=cex_line,
                    min_edge=min_edge)
    ## If there is only one point, then add a dot and label, then return directly.
    nCategory <- n
    if (inherits(n, "character")) {
        nCategory <- length(n)
    }
    if(nCategory == 1) {
        p <- ggraph(g,"tree") + geom_node_point(color="red", size=5) +
               geom_node_text(aes_(label=~name))
        return(p)
    }
    p <- adj_layout(g = g, layout = layout, coords = coords)
    
    ## add edge
    if (with_edge & length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }

    # ggData <- p$data
    # if show group cricle or group label, Process p$data and assign color to the group label
    if (group_category || node_label == "all" || node_label == "group") {         
        ggData <- groupNode(p = p, enrichDf = y, nWords = nWords, 
            clusterFunction =  clusterFunction, nCluster = nCluster)
        p$data <- ggData  
    }

    ## if group_category, add circles
    if (group_category) {
         p <- add_ellipse(p = p, group_legend = group_legend, 
            label_style = label_style, ...)
    }

    ## add dot
    p <- add_category_nodes(p = p, cex_category = cex_category, color = color)
    ## add node label
    if (node_label == "all" || node_label == "category") 
        p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext)
    ## add group label
    if (node_label == "all" || node_label == "group") {   
        label_location <- get_label_location(ggData = ggData, label_format = label_format)
        p <- add_group_label(label_style = label_style, repel = repel, shadowtext = shadowtext, p = p,
            label_location = label_location, label_group = label_group,
            cex_label_group = cex_label_group)
    }
    p + coord_equal() + 
        guides(size  = guide_legend(order = 1), 
               fill = guide_colorbar(order = 2))
}




##' @rdname emapplot
##' @importFrom igraph E "E<-"
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggplot2 coord_equal
##' @importFrom ggplot2 labs
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom scatterpie geom_scatterpie
##' @importFrom scatterpie geom_scatterpie_legend
##' @importClassesFrom DOSE compareClusterResult
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'
##' @param legend_n number of circle in legend
##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis, 
##' the default value is 1.
##' @importFrom stats setNames
emapplot.compareClusterResult <- function(x, showCategory = 30,
                                          layout = NULL,
                                          coords = NULL,
                                          split = NULL, pie = "equal",
                                          legend_n = 5, cex_category = 1,
                                          cex_line = 1, min_edge=0.2, 
                                          cex_label_category  = 1, 
                                          shadowtext = TRUE, 
                                          with_edge = TRUE,
                                          group_category = FALSE, 
                                          label_format = 30,
                                          group_legend = FALSE,
                                          node_label  = "category",
                                          label_style = "shadowtext", 
                                          repel = FALSE, cex_label_group = 1,
                                          nWords = 4, 
                                          clusterFunction = stats::kmeans,
                                          nCluster = NULL, 
                                          cex_pie2axis = 1, 
                                          ...) {
                                       
    has_pairsim(x)
    label_size_category <- 3
    label_group <- 3
    # y <- get_selected_category(showCategory, x, split)
    y <- fortify(x, showCategory = showCategory,
                 includeAll = TRUE, split = split)
    y$Cluster <- sub("\n.*", "", y$Cluster)

    if ("core_enrichment" %in% colnames(y)) { ## for GSEA result
        y$geneID <- y$core_enrichment
    }
    ## Data structure transformation, combining the same ID (Description) genes
    mergedEnrichDf <- merge_compareClusterResult(y)
     
    ## get ggraph object and add edge
    p <- build_ggraph(x = x, enrichDf = y, mergedEnrichDf = mergedEnrichDf, cex_category = cex_category, 
        pie = pie, layout = layout, coords = coords, cex_line=cex_line,
                        min_edge=min_edge, pair_sim = x@termsim,
                        method = x@method, with_edge = with_edge)
    if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(mergedEnrichDf)) | nrow(mergedEnrichDf) == 1)
        return(p)

    # ggData <- p$data
    # if show group cricle or group label, Process p$data and assign color to the group label
    if (group_category || node_label == "all" || node_label == "group") {    
        ggData <- groupNode(p = p, enrichDf = y, nWords = nWords, 
            clusterFunction =  clusterFunction, nCluster = nCluster)
        p$data <- ggData
    }      
    ## add circle
    if (group_category) {
        p <- add_ellipse(p = p, group_legend = group_legend, 
            label_style = label_style, ...)
    }
       
    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- get_pie_data(enrichDf = y, pie = pie, mergedEnrichDf = mergedEnrichDf, cex_pie2axis = cex_pie2axis, 
                                   p = p, cex_category = cex_category)


    
    ## add dot and node label
    p <- add_pie_node(p = p, ID_Cluster_mat = ID_Cluster_mat, 
                  node_label = node_label, cex_category = cex_category,
                  cex_pie2axis = cex_pie2axis, 
                  cex_label_category = cex_label_category,
                  shadowtext = shadowtext, legend_n = legend_n,
                  label_size_category = label_size_category)
    ## add group label
    if (node_label == "all" || node_label == "group") {   
        label_location <- get_label_location(ggData = ggData, label_format = label_format)
        p <- add_group_label(label_style = label_style, repel = repel, shadowtext = shadowtext, p = p,
            label_location = label_location, label_group = label_group,
            cex_label_group = cex_label_group)
    }    
    return(p)
}

