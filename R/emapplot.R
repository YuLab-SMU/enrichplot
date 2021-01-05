##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "enrichResult"),
          function(x, showCategory = 30, color = "p.adjust",
                   layout = "nicely", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, color = "p.adjust",
                   layout = "nicely", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 30, color = "p.adjust",
                   layout = "nicely", ...) {

              emapplot.compareClusterResult(x, showCategory = showCategory,
                                            color=color, layout = layout, ...)
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
##' @param cex_line scale of line width
##' @param min_edge minimum percentage of overlap genes to display the edge,
##' should between 0 and 1, default value is 0.2
##' @param cex_label_category scale of category node label size
##' @param cex_category number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @param shadowtext a logical value, whether to use shadow font.
##' @author Guangchuang Yu
emapplot.enrichResult <- function(x, showCategory = 30, color="p.adjust",
                                  layout = "nicely", min_edge=0.2,
                                  cex_label_category  = 1, cex_category = 1,
                                  cex_line = 1, shadowtext = TRUE) {
    has_pairsim(x)
    # if (!is.null(node_label_size))
    #     message("node_label_size parameter has been changed to 'cex_label_category'")
    # if (!is.null(node_scale))
    #     message("node_scale parameter has been changed to 'cex_category'")
    # if (is.null(cex_category)) {
    #     if (!is.null(node_scale)) {
    #         cex_category <- node_scale
    #     } else {
    #         cex_category <- 1
    #     }
    # }

    # if (!is.null(line_scale))
    #     message("line_scale parameter has been changed to 'cex_line'")
    # if (is.null(cex_line)) {
    #     if (!is.null(line_scale)) {
    #         cex_line <- line_scale
    #     } else {
    #         cex_line <- 1
    #     }
    # }
    label_category <- 5
    n <- update_n(x, showCategory)
    # geneSets <- geneInCategory(x) ## use core gene for gsea result


    y <- as.data.frame(x)
    g <- get_igraph(x=x, y=y, n=n, color=color, cex_line=cex_line,
                    min_edge=min_edge)
    if(n == 1) {
        return(ggraph(g,"tree") + geom_node_point(color="red", size=5) +
               geom_node_text(aes_(label=~name)))
    }

    p <- ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }
    p <- p + geom_node_point(aes_(color=~color, size=~size)) 
    p <- add_node_label(p, data = NULL, label_category,
        cex_label_category, shadowtext)
    p + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color,
                               guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8) * cex_category)

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
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @param legend_n number of circle in legend
##' @param cex_line scale of line width
##' @param min_edge minimum percentage of overlap genes to display the edge, should between 0 and 1, default value is 0.2
##' @importFrom stats setNames
emapplot.compareClusterResult <- function(x, showCategory = 30,
                                          color = "p.adjust", layout = "nicely",
                                          split = NULL, pie = "equal",
                                          legend_n = 5, cex_category = 1,
                                          cex_line = 1, min_edge=0.2, 
                                          cex_label_category  = 1, 
                                          shadowtext = TRUE) {
    has_pairsim(x)
    # if (!is.null(node_label_size))
    #     message("node_label_size parameter has been changed to 'cex_label_category'")

    # if (!is.null(pie_scale))
    #     message("pie_scale parameter has been changed to 'cex_category'")

    # if (is.null(cex_category)) {
    #     if (!is.null(pie_scale)) {
    #         cex_category <- pie_scale
    #     } else {
    #         cex_category <- 1
    #     }
    # }

    label_category <- 3
    ## pretreatment of x, just like dotplot do
    ## If showCategory is a number, keep only the first showCategory of each group
    ## Otherwise keep the total showCategory rows
    y <- get_selected_category(showCategory, x, split)
    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- merge_compareClusterResult(y)

    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/",
                                  fixed = TRUE), y_union$ID)

    g <- build_emap_graph(y=y_union,geneSets=geneSets,color=color,
                          cex_line=cex_line, min_edge=min_edge,
                          pair_sim = x@termsim, method = x@method)

    p <- build_ggraph(y = y, g = g, y_union = y_union, cex_category = cex_category,
               pie = pie, layout = layout)               
    if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(y_union)) | nrow(y_union) == 1)
        return(p)
    p <- ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }

    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- prepare_pie_category(y,pie=pie)

    # plot the edge
    # get the X-coordinate and y-coordinate of pies
    aa <- p$data
    desc <- y_union$Description[match(rownames(ID_Cluster_mat),
                                      y_union$Description)]
    i <- match(desc, aa$name)
    ID_Cluster_mat$x <- aa$x[i]
    ID_Cluster_mat$y <- aa$y[i]

    #Change the radius value to fit the pie plot
    radius <- NULL
    ID_Cluster_mat$radius <- sqrt(aa$size[i] / sum(aa$size) * cex_category)
    #ID_Cluster_mat$radius <- sqrt(aa$size / pi)

    x_loc1 <- min(ID_Cluster_mat$x)
    y_loc1 <- min(ID_Cluster_mat$y)
    ## x_loc2 <- min(ID_Cluster_mat$x)
    ## y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))
    if(ncol(ID_Cluster_mat) > 4) {
        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
            cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
            coord_equal() 
            # geom_node_text(aes_(label=~name), repel=TRUE,
            # size = label_category * cex_label_category, bg.color = "white") + 
        p <- add_node_label(p, data = NULL, label_category,
            cex_label_category, shadowtext)
        p <- p + theme_void() +
            geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
                n = legend_n,
                labeller=function(x) round(sum(aa$size) * x^2 / cex_category)) +
            labs(fill = "Cluster")
        return(p)
    }
    ## annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")
    title <- colnames(ID_Cluster_mat)[1]
    # p + geom_node_point(aes_(color=~color, size=~size))
    # p + geom_node_text(aes_(label=~name), repel=TRUE,
    #     size = label_category * cex_label_category, bg.color = "white") + 
    p <- add_node_label(p, data = NULL, label_category,
            cex_label_category, shadowtext)
    p <- p + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color,
                               guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8) * cex_category)  +labs(title= title)
}

