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
##' 'randomly', 'fr', 'kk', 'drl' or 'lgl'. Will be removed in the next version.
##' Will be removed in the next version.
##' @param coords a data.frame with two columns: 'x' for X-axis coordinate and 'y' for Y-axis coordinate.
##' Will be removed in the next version.
##' @param color Variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @param cex_line Scale of line width.
##' Will be removed in the next version.
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##' Will be removed in the next version.
##' @param cex_label_category Scale of category node label size.
##' Will be removed in the next version.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' Will be removed in the next version.
##' @param shadowtext a logical value, whether to use shadow font.
##' @param label_style style of group label, one of "shadowtext" and "ggforce".
##' Will be removed in the next version.
##' @param repel whether to correct the position of the label. Defaults to FALSE.
##' @param node_label Select which labels to be displayed,
##' one of 'category', 'group', 'all' and 'none'.
##' @param with_edge Logical, if TRUE, draw the edges of the network diagram.
##' Will be removed in the next version.
##' @param group_category a logical, if TRUE, group the category.
##' Will be removed in the next version.
##' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
##' The default is FALSE.
##' Will be removed in the next version.
##' @param cex_label_group Numeric, scale of group labels size, the default value is 1.
##' Will be removed in the next version.
##' @param nWords Numeric, the number of words in the cluster tags, the default value is 4.
##' Will be removed in the next version.
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' Will be removed in the next version.
##' @param clusterFunction function of Clustering method, such as stats::kmeans(the default),
##' cluster::clara, cluster::fanny or cluster::pam.
##' Will be removed in the next version.
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' Will be removed in the next version.
##' @param layout.params list, the parameters to control the layout.
##' see the layout.params in the following.
##' layout.params control the attributes of layout, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{layout} Layout of the map, e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
##'         'randomly', 'fr', 'kk', 'drl' or 'lgl'..
##'         \item \code{coords} a data.frame with two columns: 'x' for X-axis coordinate and 
##'         'y' for Y-axis coordinate.
##'     }
##' @param edge.params list, the parameters to control the edge.
##' see the edge.params in the following.
##' edge.params control the attributes of edge, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{show} Logical, if TRUE (the default), draw the edges of the network diagram.
##'         \item \code{min} The minimum similarity threshold for whether 
##'         two nodes are connected, should between 0 and 1, default value is 0.2.
##'     }
##' @param cex.params list, the parameters to control the edge.
##' see the cex.params in the following.
##' cex.params control the attributes of edge, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{category_node} Number indicating the amount by which plotting category
##'         nodes should be scaled relative to the default.
##'         \item \code{category_label} Scale of category node label size. 
##'         \item \code{line} Scale of line width.
##'         \item \code{pie2axis} It is used to adjust the relative size of the pie chart on the coordinate axis, 
##'         the default value is 1.
##'         \item \code{label_group} Numeric, scale of group labels size, the default value is 1.
##'     }
##' @param cluster.params list, the parameters to control the attributes of highlighted nodes and edges.
##' see the cluster.params in the following.
##' cluster.params control the attributes of highlight, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{cluster} a logical, if TRUE, group the category.
##'         \item \code{method} function of Clustering method, such as stats::kmeans(the default),
##'         cluster::clara, cluster::fanny or cluster::pam.
##'         \item \code{n} Numeric, the number of clusters, 
##'         the default value is square root of the number of nodes.
##'         \item \code{legend} Logical, if TRUE, the grouping legend will be displayed.
##'         The default is FALSE.
##'         \item \code{label_style} style of group label, one of "shadowtext" and "ggforce".
##'         \item \code{label_words_n} Numeric, the number of words in the cluster tags, the default value is 4.
##'         \item \code{label_format} a numeric value sets wrap length, alternatively a
##'         custom function to format axis labels.
##'     }
##' @param hilight.params list, the parameters to control the attributes of highlighted nodes and edges.
##' see the hilight.params in the following.
##' hilight.params control the attributes of highlight, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{category} category nodes to be highlight.
##'         \item \code{alpha_hilight} alpha of highlighted nodes.
##'         \item \code{alpha_no_hilight} alpha of unhighlighted nodes.
##'     }
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
emapplot.enrichResult <- function(x, 
                                  showCategory = 30, 
                                  layout,                    # removed
                                  coords,                    # removed
                                  color = "p.adjust",     
                                  min_edge,                  # removed
                                  cex_label_category,        # removed
                                  cex_category,              # removed
                                  cex_line,                  # removed
                                  shadowtext = TRUE,      
                                  label_style,               # removed
                                  repel,                     # removed       
                                  node_label  = "category",     
                                  with_edge,                 # removed
                                  group_category,            # removed
                                  group_legend,              # removed                      
                                  cex_label_group,           # removed
                                  nWords,                    # removed
                                  label_format,              # removed
                                  clusterFunction,           # removed
                                  nCluster,                  # removed
                                  layout.params = list(
                                      layout = NULL, 
                                      coords = NULL             
                                  ),

                                  edge.params = list(
                                      show = TRUE,              
                                      min = 0.2                           
                                  ),
                                  cex.params = list(
                                      category_node = 1,        
                                      category_label = 1,       
                                      line = 1                  
                                  ),
                                  hilight.params = list(
                                      category = NULL,          
                                      alpha_hilight = 1,        
                                      alpha_no_hilight = 0.3    
                                  ),  
                                  cluster.params = list(
                                      cluster = FALSE,          
                                      method = stats::kmeans,   
                                      n = NULL,                               
                                      legend = FALSE,           
                                      label_style = "shadowtext",
                                      label_words_n = 4,        
                                      label_format = 30         
                                  ),
                                  ...) {
    has_pairsim(x)
    label_size_category <- 5
    label_group <- 3

    # change parameter name
    ##############################################################
    params_df <- as.data.frame(rbind(
        c("layout", "layout.params", "layout"),
        c("coords", "layout.params", "coords"),

        c("with_edge", "edge.params", "show"),
        c("min_edge", "edge.params", "min"),

        c("cex_category", "cex.params", "category_node"),
        c("cex_label_category", "cex.params", "category_label"),
        c("cex_line", "cex.params", "line"),

        c("group_category", "cluster.params", "cluster"),
        c("clusterFunction", "cluster.params", "method"),
        c("nCluster", "cluster.params", "n"),
        c("group_legend", "cluster.params", "legend"),
        c("label_style", "cluster.params", "label_style"),
        c("nWords", "cluster.params", "label_words_n"),
        c("label_format", "cluster.params", "label_format"))
    )
    colnames(params_df) <- c("original", "listname", "present")
    rownames(params_df) <- params_df$original

 
    # layout.params
    default.layout.params <- list(
        layout = NULL,                       
        coords = NULL            
    )
    layout.params <- reset_params(defaultp=default.layout.params, 
                                   inputp=enquo(layout.params))

                           
    if (!missing(layout)) {
        warn <- get_param_change_message("layout", params_df)
        warning(warn)
        layout.params[[params_df["layout", "present"]]] <- layout
    }  

    if (!missing(coords)) {
        warn <- get_param_change_message("coords", params_df)
        warning(warn)
        layout.params[[params_df["coords", "present"]]] <- coords
    }  
    layout <- layout.params[["layout"]]
    coords <- layout.params[["coords"]]


    # edge.params
    default.edge.params <- list(
        show = TRUE,                  
        min = 0.2                             
    )
    edge.params <- reset_params(defaultp=default.edge.params, 
                                   inputp=enquo(edge.params))

    if (!missing(with_edge)) {
        warn <- get_param_change_message("with_edge", params_df)
        warning(warn)
        edge.params[[params_df["with_edge", "present"]]] <- with_edge
    }  

    if (!missing(min_edge)) {
        warn <- get_param_change_message("min_edge", params_df)
        warning(warn)
        edge.params[[params_df["min_edge", "present"]]] <- min_edge
    }  

    with_edge <- edge.params[["show"]]
    min_edge <- edge.params[["min"]]

    # cex.params
    default.cex.params <- list(
        category_node = 1,          
        category_label = 1,         
        line = 1                                       
    )
    cex.params <- reset_params(defaultp=default.cex.params, 
                                   inputp=enquo(cex.params))

    if (!missing(cex_category)) {
        warn <- get_param_change_message("cex_category", params_df)
        warning(warn)
        cex.params[[params_df["cex_category", "present"]]] <- cex_category
    }  

    if (!missing(cex_label_category)) {
        warn <- get_param_change_message("cex_label_category", params_df)
        warning(warn)
        cex.params[[params_df["cex_label_category", "present"]]] <- cex_label_category
    }

    if (!missing(cex_line)) {
        warn <- get_param_change_message("cex_line", params_df)
        warning(warn)
        cex.params[[params_df["cex_line", "present"]]] <- cex_line
    }   
    cex_category <- cex.params[["category_node"]]
    cex_label_category <- cex.params[["category_label"]]
    cex_line <- cex.params[["line"]]

    # cluster.params
    default.cluster.params <- list(
        cluster = FALSE,            
        method = stats::kmeans,     
        n = NULL,                             
        legend = FALSE,             
        label_style = "shadowtext", 
        label_words_n = 4,          
        label_format = 30                          
    )
    cluster.params <- reset_params(defaultp=default.cluster.params, 
                                   inputp=enquo(cluster.params))

    if (!missing(group_category)) {
        warn <- get_param_change_message("group_category", params_df)
        warning(warn)
        cluster.params[[params_df["group_category", "present"]]] <- group_category
    } 

    if (!missing(clusterFunction)) {
        warn <- get_param_change_message("clusterFunction", params_df)
        warning(warn)
        cluster.params[[params_df["clusterFunction", "present"]]] <- clusterFunction
    }  

    if (!missing(nCluster)) {
        warn <- get_param_change_message("nCluster", params_df)
        warning(warn)
        cluster.params[[params_df["nCluster", "present"]]] <- nCluster
    }  

    if (!missing(group_legend)) {
        warn <- get_param_change_message("group_legend", params_df)
        warning(warn)
        cluster.params[[params_df["group_legend", "present"]]] <- group_legend
    }  

    if (!missing(label_style)) {
        warn <- get_param_change_message("label_style", params_df)
        warning(warn)
        cluster.params[[params_df["label_style", "present"]]] <- label_style
    }  

    if (!missing(nWords)) {
        warn <- get_param_change_message("nWords", params_df)
        warning(warn)
        cluster.params[[params_df["nWords", "present"]]] <- nWords
    }  

    if (!missing(label_format)) {
        warn <- get_param_change_message("label_format", params_df)
        warning(warn)
        cluster.params[[params_df["label_format", "present"]]] <- label_format
    }    

    group_category <- cluster.params[["cluster"]]
    clusterFunction <- cluster.params[["method"]]
    nCluster <- cluster.params[["n"]]
    group_legend <- cluster.params[["legend"]]
    label_style <- cluster.params[["label_style"]]
    nWords <- cluster.params[["label_words_n"]]
    label_format <- cluster.params[["label_format"]]

    ##########################################################

    default.hilight.params <- list(
        category = NULL,
        alpha_hilight = 1,
        alpha_no_hilight = 0.3
    )
    hilight.params <- reset_params(defaultp=default.hilight.params, 
                                   inputp=enquo(hilight.params))
    hilight_category <- hilight.params[["category"]]
    alpha_hilight <- hilight.params[["alpha_hilight"]]
    alpha_nohilight <- hilight.params[["alpha_no_hilight"]]
    n <- update_n(x, showCategory)
    y <- as.data.frame(x)
    ## get graph.data.frame() object
    g <- get_igraph(x=x, nCategory=n, color=color, cex_line=cex_line,
                    min_edge=min_edge)
    hilight_category <- intersect(hilight_category, attr(V(g), "names"))
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
    # if (!is.null(hilight_category) && length(hilight_category) > 0) {
    #     edges <- attr(E(g), "vnames")
    #     E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
    #     hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
    #     E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
    #     # E(g)$alpha[hilight_edge] <- alpha_hilight
    # } else {
    #     E(g)$alpha <- rep(min(0.8, alpha_hilight), length(E(g)))
    # }
    g <- edge_add_alpha(g, hilight_category, alpha_nohilight, alpha_hilight)
    p <- adj_layout(g = g, layout = layout, coords = coords)
    p$data$alpha <- rep(1, nrow(p$data))
    if (!is.null(hilight_category) && length(hilight_category) > 0) {
        if (length(hilight_category) > 0) {
            p$data$alpha <- rep(alpha_nohilight, nrow(p$data))
            ii <- match(hilight_category, p$data$name)
            p$data$alpha[ii] <- alpha_hilight
        }
    }
    ## add edge
   
    if (with_edge & length(E(g)$width) > 0) {
        p <- p + geom_edge_link(aes_(width=~I(width), alpha=~I(alpha)),
                                colour='darkgrey', show.legend = FALSE)
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
               fill = guide_colorbar(order = 2),
               alpha = "none")
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
##' Will be removed in the next version.
##' @param legend_n number of circle in legend
##' Will be removed in the next version.
##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis, 
##' the default value is 1.
##' Will be removed in the next version.
##' @param pie.params list, the parameters to control the attributes of pie nodes.
##' see the pie.params in the following.
##' pie.params control the attributes of pie nodes, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{pie} proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'.
##'         \item \code{legend_n} number of circle in legend.
##'     }
##' @importFrom stats setNames
emapplot.compareClusterResult <- function(x, 
                                          showCategory = 30,
                                          layout,                   # removed
                                          coords,                   # removed
                                          split = NULL,      
                                          pie,                   # removed
                                          legend_n,                    # removed
                                          cex_category,                # removed
                                          cex_line,                    # removed
                                          min_edge,                    # removed
                                          cex_label_category,         # removed
                                          shadowtext,      
                                          with_edge,                # removed
                                          group_category,          # removed
                                          label_format,               # removed
                                          group_legend,            # removed
                                          node_label  = "category",      
                                          label_style,      # removed
                                          repel = FALSE,      
                                          cex_label_group,             # removed
                                          nWords,                      # removed
                                          clusterFunction, # removed
                                          nCluster,                 # removed
                                          cex_pie2axis,                # removed

                                          pie.params = list(
                                              pie = "equal",             
                                              legend_n = 5               
                                          ),    
                                          layout.params = list(
                                              layout = NULL,             
                                              coords = NULL              
                                          ),
                                          edge.params = list(
                                              show = TRUE,               
                                              min = 0.2                  
                                          ),
                                          cluster.params = list(
                                              cluster = FALSE,           
                                              method = stats::kmeans,    
                                              n = NULL,                  
                                              legend = FALSE,            
                                              label_style = "shadowtext",
                                              label_words_n = 4,         
                                              label_format = 30          
                                          ),
                                          cex.params = list(
                                              category_node = 1,        
                                              category_label = 1,       
                                              line = 1,                 
                                              pie2axis = 1,             
                                              label_group = 1         
                                          ),
                                          hilight.params=list(
                                              category = NULL,          
                                              alpha_hilight = 1,          
                                              alpha_no_hilight = 0.3      
                                          ),  
                                          ...) {
                                       
    has_pairsim(x)
    label_size_category <- 3
    label_group <- 3

    # change parameter name
    ##############################################################
    params_df <- as.data.frame(rbind(
        c("pie", "pie.params", "pie"),
        c("legend_n", "pie.params", "legend_n"),

        c("layout", "layout.params", "layout"),
        c("coords", "layout.params", "coords"),

        c("min_edge", "edge.params", "show"),
        c("with_edge", "edge.params", "min"),

        c("group_category", "cluster.params", "cluster"),
        c("clusterFunction", "cluster.params", "method"),
        c("nCluster", "cluster.params", "n"),
        c("label_style", "cluster.params", "legend"),
        c("group_legend", "cluster.params", "label_style"),
        c("nWords", "cluster.params", "label_words_n"),
        c("label_format", "cluster.params", "label_format"),

        c("cex_category", "cex.params", "category_node"),
        c("cex_label_category", "cex.params", "category_label"),
        c("cex_line", "cex.params", "line"),
        c("cex_pie2axis", "cex.params", "pie2axis"),
        c("cex_label_group", "cex.params", "label_group"))
    )
    colnames(params_df) <- c("original", "listname", "present")
    rownames(params_df) <- params_df$original

 
    # pie.params
    default.pie.params <- list(
        pie = "equal",          
        legend_n = 5                
    )
    pie.params <- reset_params(defaultp=default.pie.params, 
                               inputp=enquo(pie.params))

                           
    if (!missing(pie)) {
        warn <- get_param_change_message("pie", params_df)
        warning(warn)
        pie.params[[params_df["pie", "present"]]] <- pie
    }  

    if (!missing(legend_n)) {
        warn <- get_param_change_message("legend_n", params_df)
        warning(warn)
        pie.params[[params_df["legend_n", "present"]]] <- legend_n
    }  

    pie <- pie.params[["pie"]]
    legend_n <- pie.params[["legend_n"]]


    # layout.params
    default.layout.params <- list(
        layout = NULL,                       
        coords = NULL                  
    )
    layout.params <- reset_params(defaultp=default.layout.params, 
                               inputp=enquo(layout.params))

                           
    if (!missing(layout)) {
        warn <- get_param_change_message("layout", params_df)
        warning(warn)
        layout.params[[params_df["layout", "present"]]] <- layout
    }  

    if (!missing(coords)) {
        warn <- get_param_change_message("coords", params_df)
        warning(warn)
        layout.params[[params_df["coords", "present"]]] <- coords
    }  

    layout <- layout.params[["layout"]]
    coords <- layout.params[["coords"]]

    # edge.params
    default.edge.params <- list(
        show = TRUE,                 
        min = 0.2                       
    )
    edge.params <- reset_params(defaultp=default.edge.params, 
                               inputp=enquo(edge.params))

                           
    if (!missing(with_edge)) {
        warn <- get_param_change_message("with_edge", params_df)
        warning(warn)
        edge.params[[params_df["with_edge", "present"]]] <- with_edge
    }  

    if (!missing(min_edge)) {
        warn <- get_param_change_message("min_edge", params_df)
        warning(warn)
        edge.params[[params_df["min_edge", "present"]]] <- min_edge
    }  

    with_edge <- edge.params[["show"]]
    min_edge <- edge.params[["min"]]

    # cluster.params
    default.cluster.params <- list(
        cluster = FALSE,          
        method = stats::kmeans,   
        n = NULL,                  
        legend = FALSE,           
        label_style = "shadowtext",
        label_words_n = 4,        
        label_format = 30                     
    )
    cluster.params <- reset_params(defaultp=default.cluster.params, 
                               inputp=enquo(cluster.params))

                           
    if (!missing(group_category)) {
        warn <- get_param_change_message("group_category", params_df)
        warning(warn)
        cluster.params[[params_df["group_category", "present"]]] <- group_category
    }  

    if (!missing(clusterFunction)) {
        warn <- get_param_change_message("clusterFunction", params_df)
        warning(warn)
        cluster.params[[params_df["clusterFunction", "present"]]] <- clusterFunction
    }  

    if (!missing(nCluster)) {
        warn <- get_param_change_message("nCluster", params_df)
        warning(warn)
        cluster.params[[params_df["nCluster", "present"]]] <- nCluster
    }  

    if (!missing(label_style)) {
        warn <- get_param_change_message("label_style", params_df)
        warning(warn)
        cluster.params[[params_df["label_style", "present"]]] <- label_style
    }  

    if (!missing(group_legend)) {
        warn <- get_param_change_message("group_legend", params_df)
        warning(warn)
        cluster.params[[params_df["group_legend", "present"]]] <- group_legend
    }  

    if (!missing(nWords)) {
        warn <- get_param_change_message("nWords", params_df)
        warning(warn)
        cluster.params[[params_df["nWords", "present"]]] <- nWords
    }  

    if (!missing(label_format)) {
        warn <- get_param_change_message("label_format", params_df)
        warning(warn)
        cluster.params[[params_df["label_format", "present"]]] <- label_format
    }  
    group_category <- cluster.params[["cluster"]]
    clusterFunction <- cluster.params[["method"]]
    nCluster <- cluster.params[["n"]]
    label_style <- cluster.params[["legend"]]
    group_legend <- cluster.params[["label_style"]]
    nWords <- cluster.params[["label_words_n"]]
    label_format <- cluster.params[["label_format"]]

    # cex.params
    default.cex.params <- list(
        category_node = 1,                 
        category_label = 1,           
        line = 1,                     
        pie2axis = 1,                 
        label_group = 1                        
    )
    cex.params <- reset_params(defaultp=default.cex.params, 
                               inputp=enquo(cex.params))

                           
    if (!missing(cex_category)) {
        warn <- get_param_change_message("cex_category", params_df)
        warning(warn)
        cex.params[[params_df["cex_category", "present"]]] <- cex_category
    }  

    if (!missing(cex_label_category)) {
        warn <- get_param_change_message("cex_label_category", params_df)
        warning(warn)
        cex.params[[params_df["cex_label_category", "present"]]] <- cex_label_category
    }  

    if (!missing(cex_line)) {
        warn <- get_param_change_message("cex_line", params_df)
        warning(warn)
        cex.params[[params_df["cex_line", "present"]]] <- cex_line
    }  

    if (!missing(cex_pie2axis)) {
        warn <- get_param_change_message("cex_pie2axis", params_df)
        warning(warn)
        cex.params[[params_df["cex_pie2axis", "present"]]] <- cex_pie2axis
    }  

    if (!missing(cex_label_group)) {
        warn <- get_param_change_message("cex_label_group", params_df)
        warning(warn)
        cex.params[[params_df["cex_label_group", "present"]]] <- cex_label_group
    }  


    cex_category <- cex.params[["category_node"]]
    cex_label_category <- cex.params[["category_label"]]
    cex_line <- cex.params[["line"]]
    cex_pie2axis <- cex.params[["pie2axis"]]
    cex_label_group <- cex.params[["label_group"]]



    default.hilight.params <- list(
        category = NULL,
        alpha_hilight = 1,
        alpha_no_hilight = 0.3
    )
    hilight.params <- reset_params(defaultp=default.hilight.params, 
                                   inputp=enquo(hilight.params))
    hilight_category <- hilight.params[["category"]]
    alpha_hilight <- hilight.params[["alpha_hilight"]]
    alpha_nohilight <- hilight.params[["alpha_no_hilight"]]
    # y <- get_selected_category(showCategory, x, split)
    y <- fortify(x, showCategory = showCategory,
                 includeAll = TRUE, split = split)
    y$Cluster <- sub("\n.*", "", y$Cluster)

    if ("core_enrichment" %in% colnames(y)) { ## for GSEA result
        y$geneID <- y$core_enrichment
    }
    hilight_category <- intersect(hilight_category, y$Description)
    ## Data structure transformation, combining the same ID (Description) genes
    mergedEnrichDf <- merge_compareClusterResult(y)
     
    ## get ggraph object and add edge
    p <- build_ggraph(x = x, enrichDf = y, mergedEnrichDf = mergedEnrichDf, cex_category = cex_category, 
        pie = pie, layout = layout, coords = coords, cex_line=cex_line,
                        min_edge=min_edge, pair_sim = x@termsim,
                        method = x@method, with_edge = with_edge,
                        hilight_category = hilight_category,
                        alpha_hilight = alpha_hilight,
                        alpha_nohilight = alpha_nohilight)
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
    return(p + guides(alpha = "none"))
}

