##' @rdname treeplot
##' @exportMethod treeplot
setMethod("treeplot", signature(x = "enrichResult"),
    function(x, ...) {
        treeplot.enrichResult(x, ...)
    })

##' @rdname treeplot
##' @exportMethod treeplot
setMethod("treeplot", signature(x = "gseaResult"),
    function(x, ...) {
        treeplot.enrichResult(x, ...)
    })

##' @rdname treeplot
##' @exportMethod treeplot
setMethod("treeplot", signature(x = "compareClusterResult"),
    function(x, ...) {
        treeplot.compareClusterResult(x, ...)
    })



##' @rdname treeplot
##' @param nWords The number of words in the cluster tags.
##' Will be removed in the next version.
##' @param nCluster The number of clusters, the default value is 5.
##' Will be removed in the next version.
##' @param cex_category Number indicating the amount by which plotting category.
##' nodes should be scaled relative to the default.
##' Will be removed in the next version.
##' @param label_format_cladelab label_format for group labels, a numeric value sets wrap length, 
##' alternatively a custom function to format axis labels.
##' Will be removed in the next version.
##' @param label_format_tiplab label_format for tiplabs, a numeric value sets wrap length, 
##' alternatively a custom function to format axis labels.
##' Will be removed in the next version.
##' @param offset_tiplab tiplab offset, rel object or numeric value, the bigger the number, 
##' the farther the distance between the node and the branch. 
##' The default is rel(1), when geneClusterPanel = "pie", meaning 1 *  max_radius_of_the_pies; 
##' when geneClusterPanel = "heatMap", meaning 1 * 0.16 * column_number_of_heatMap * x_range_of_tree;
##' when geneClusterPanel = "dotplot", meaning 1 * 0.09 * column_number_of_dotplot * x_range_of_tree.
##' Will be removed in the next version.
##' @param offset rel object or numeric value, distance bar and tree,
##' offset of bar and text from the clade, default is rel(1),
##' meaning 1 * 1.2 * x_range_of_tree plus distance_between_tree_and_tiplab
##' (1 * (1.2 * x_range_of_tree + distance_between_tree_and_tiplab)).
##' Will be removed in the next version.
##' @param fontsize The size of text, default is 4.
##' @param hclust_method Method of hclust. This should be (an unambiguous abbreviation of) one of "ward.D", 
##' "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
##' Will be removed in the next version.
##' @param group_color A vector of group colors, the length of the vector should be the same as nCluster.
##' Will be removed in the next version.
##' @param extend Numeric, extend the length of bar, default is 0.3.
##' Will be removed in the next version.
##' @param hilight Logical value, if TRUE(default), add ggtree::geom_hilight() layer.
##' Will be removed in the next version.
##' @param hexpand expand x limits by amount of xrange * hexpand.
##' Will be removed in the next version.
##' @param align control the align direction of the edge of high light rectangular.
##' Options is 'none', 'left', 'right', 'both (default)'.
##' Will be removed in the next version.
##' @param hilight.params list, the parameters to control the attributes of highlight layer.
##' see the hilight.params in the following.
##' hilight.params control the attributes of highlight layer, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{hilight} Logical value, if TRUE(default), add ggtree::geom_hilight() layer.
##'         \item \code{align} control the align direction of the edge of high light rectangular.
##'         Options is 'none', 'left', 'right', 'both (default)'.
##'     }
##' @param offset.params list, the parameters to control the offset.
##' see the offset.params in the following.
##' offset.params control the attributes of offset, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{bar_tree} rel object or numeric value, distance bar and tree,
##'         offset of bar and text from the clade, default is rel(1),
##'         meaning 1 * 1.2 * x_range_of_tree plus distance_between_tree_and_tiplab
##'         (1 * (1.2 * x_range_of_tree + distance_between_tree_and_tiplab)).
##'         \item \code{tiplab} tiplab offset, rel object or numeric value, the bigger the number, 
##'         the farther the distance between the node and the branch. 
##'         The default is rel(1), when clusterPanel = "pie", meaning 1 *  max_radius_of_the_pies; 
##'         when clusterPanel = "heatMap", meaning 1 * 0.16 * column_number_of_heatMap * x_range_of_tree;
##'         when clusterPanel = "dotplot", meaning 1 * 0.09 * column_number_of_dotplot * x_range_of_tree.
##'         \item \code{extend} Numeric, extend the length of bar, default is 0.3.
##'         \item \code{hexpand} expand x limits by amount of xrange * hexpand.
##'     }
##' @param cluster.params list, the parameters to control the attributes of highlighted nodes and edges.
##' see the cluster.params in the following.
##' cluster.params control the attributes of highlight, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{method} function of Clustering method, such as stats::kmeans(the default),
##'         cluster::clara, cluster::fanny or cluster::pam.
##'         \item \code{n} Numeric, the number of clusters, 
##'         the default value is square root of the number of nodes.
##'         \item \code{color} A vector of group colors, the length of the vector should be the same as nCluster.
##'         \item \code{label_words_n} Numeric, the number of words in the cluster tags, the default value is 4.
##'         \item \code{label_format} A numeric value sets wrap length, alternatively a
##'         custom function to format axis labels.
##'     }
##' @importFrom ggtree `%<+%`
##' @importFrom ggtree ggtree
##' @importFrom ggtree geom_tiplab
##' @importFrom ggtree geom_tippoint
##' @importFrom ggtree groupClade
##' @importFrom ggtree geom_cladelab
##' @importFrom ggplot2 coord_cartesian
##' @importFrom ggplot2 scale_colour_continuous
##'
treeplot.enrichResult <- function(x, 
                                  showCategory = 30,
                                  color = "p.adjust",
                                  nWords = 4,                   # removed
                                  nCluster = 5,                 # removed
                                  cex_category = 1,
                                  label_format = NULL, 
                                  label_format_cladelab = 30,   # removed
                                  label_format_tiplab = NULL, 
                                  fontsize = 4, 
                                  offset = rel(1),              # removed
                                  offset_tiplab = rel(1),       # removed
                                  hclust_method = "ward.D",     # removed
                                  group_color = NULL,           # removed
                                  extend = 0.3,                 # removed
                                  hilight = TRUE,               # removed
                                  hexpand = .1,                 # removed
                                  align = "both",               # removed
                                  hilight.params = list(
                                      hilight = TRUE,
                                      align = "both"
                                  ),
                                  offset.params = list(
                                      bar_tree = rel(1),
                                      tiplab = rel(1),
                                      extend = 0.3,
                                      hexpand = .1
                                  ),
                                  cluster.params = list(
                                      method = "ward.D",
                                      n = 5,
                                      color = NULL,
                                      label_words_n = 4,
                                      label_format = 30
                                  ),
                                  ...) {


    # change parameter name
    ##############################################################
    params_df <- as.data.frame(rbind(
        c("hilight", "hilight.params", "hilight"),
        c("align", "hilight.params", "align"),

        c("offset", "offset.params", "bar_tree"),
        c("offset_tiplab", "offset.params", "tiplab"),
        c("extend", "offset.params", "extend"),
        c("hexpand", "offset.params", "hexpand"),

        c("hclust_method", "cluster.params", "method"),
        c("nCluster", "cluster.params", "n"),
        c("group_color", "cluster.params", "color"),
        c("nWords", "cluster.params", "label_words_n"),
        c("label_format_cladelab", "cluster.params", "label_format"))
    )
    colnames(params_df) <- c("original", "listname", "present")
    rownames(params_df) <- params_df$original

    default.hilight.params <- list(
        hilight = TRUE,                        
        align = "both" 
    )
    default.offset.params <- list(
        bar_tree = rel(1),       
        tiplab = rel(1),
        extend = 0.3,          
        hexpand = .1           
    )
    default.cluster.params <- list(
        method = "ward.D",   
        n = 5,               
        color = NULL,        
        label_words_n = 4,   
        label_format = 30    
    )
    # use modifyList to change the values of parameter 
    hilight.params <- modifyList(default.hilight.params, hilight.params)
    offset.params <- modifyList(default.offset.params, offset.params)
    cluster.params <- modifyList(default.cluster.params, cluster.params)
    params_list <- list(x = x,
        showCategory = showCategory,
        color = color,
        nWords = nWords,                   
        nCluster = nCluster,                 
        cex_category = cex_category,
        label_format = label_format, 
        label_format_cladelab = label_format_cladelab,   
        label_format_tiplab = label_format_tiplab, 
        fontsize = fontsize, 
        offset = offset,              
        offset_tiplab = offset_tiplab,       
        hclust_method = hclust_method,     
        group_color = group_color,           
        extend = extend,                 
        hilight = hilight,               
        hexpand = hexpand,                 
        align = align,               
        hilight.params = hilight.params,
        offset.params = offset.params,
        cluster.params = cluster.params
    )

    # get all parameters value
    args <- as.list(match.call())
    removed_params <- intersect(params_df$original, names(args))
    if (length(removed_params) > 0) {
        for (i in removed_params) {
            params_list[[params_df[i, 2]]][[params_df[i, 3]]] <- get(i)
            warn <- get_param_change_message(i, params_df)
            warning(warn)
        }
    }

    hilight.params <- params_list[["hilight.params"]]
    offset.params <- params_list[["offset.params"]]
    cluster.params <- params_list[["cluster.params"]]


    hilight <- hilight.params[["hilight"]]
    align <- hilight.params[["align"]]
    offset <- offset.params[["bar_tree"]]
    offset_tiplab <- offset.params[["tiplab"]]
    extend <- offset.params[["extend"]]
    hexpand <- offset.params[["hexpand"]]
    hclust_method <- cluster.params[["method"]]
    nCluster <- cluster.params[["n"]]
    group_color <- cluster.params[["color"]]
    nWords <- cluster.params[["label_words_n"]]
    label_format_cladelab <- cluster.params[["label_format"]]    


    group <- p.adjust <- count<- NULL
    # to compatible with older versions
    if (!is.null(label_format)) {
        label_format_cladelab <- label_format
    }
    # if (class(x) == "gseaResult")
    if (inherits(x, "gseaResult")) {
        x@result$Count <- x$core_enrichment %>%
            strsplit(split = "/")  %>%
            vapply(length, FUN.VALUE = 1)   
    } 

    n <- update_n(x, showCategory)
    if (is.numeric(n)) {
        keep <- seq_len(n)
    } else {
        keep <- match(n, rownames(x@termsim))
    }

    if (length(keep) == 0) {
        stop("no enriched term found...")
    }
    ## Fill the upper triangular matrix completely
    termsim2 <- fill_termsim(x, keep)

    ## Use the ward.D method to avoid overlapping ancestor nodes of each group
    hc <- stats::hclust(stats::as.dist(1- termsim2),
                        method = hclust_method)
    clus <- stats::cutree(hc, nCluster)
    d <- data.frame(label = names(clus),
        #node = seq_len(length(clus)),
        color = x[keep, as.character(color)],
        count = x$Count[keep])
    
    ## Group the nodes.
    p <- group_tree(hc = hc, clus = clus, d = d, offset_tiplab = offset_tiplab, 
        nWords = nWords, label_format_cladelab = label_format_cladelab, 
        label_format_tiplab = label_format_tiplab, offset = offset, 
        fontsize = fontsize, group_color = group_color, extend = extend, 
        hilight = hilight, cex_category = cex_category, align = align, align_tiplab = FALSE,
        color = color)     
    # xlim <-  c(0, xlim * 3 * max(p$data$x))
    # p + coord_cartesian(xlim = xlim) +
    #   p + ggnewscale::new_scale_colour() +
    #     geom_tippoint(aes(color = color, size = count)) +
      p + scale_size_continuous(name = "number of genes",
                              range = c(3, 8) * cex_category) + 
        ggtree::hexpand(ratio = hexpand) + 
        guides(size  = guide_legend(order = 1), 
               color = guide_colorbar(order = 2))
}


##' @rdname treeplot
##' @param split Separate result by 'category' variable.
##' @param legend_n Number of circle in legend, the default value is 3.
##' Will be removed in the next version.
##' @param geneClusterPanel one of "heatMap"(default), "dotplot", "pie".
##' Will be removed in the next version.
##' @param pie Used only when geneClusterPanel = "pie", 
##' proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'.
##' Will be removed in the next version.
##' @param clusterPanel.params list, the parameters to control the attributes of cluster panel.
##' see the clusterPanel.params in the following.
##' clusterPanel.params control the attributes of cluster panel, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{clusterPanel} one of "heatMap"(default), "dotplot", "pie".
##'         \item \code{pie} pUsed only when ClusterPanel = "pie", 
##'         proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'.
##'         \item \code{legend_n} number of circle in legend.
##'         \item \code{colnames_angle} set the angle of colnames.
##'     }
##' @importFrom ggtree nodepie
##' @importFrom ggtree geom_inset
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom rlang check_installed
treeplot.compareClusterResult <-  function(x, 
                                      showCategory = 5,
                                      color = "p.adjust",
                                      nWords = 4,                     # removed
                                      nCluster = 5,                   # removed
                                      cex_category = 1, 
                                      split = NULL,
                                      label_format = NULL, 
                                      label_format_cladelab = 30,     # removed
                                      label_format_tiplab = NULL, 
                                      fontsize = 4, 
                                      offset = rel(1),                 # removed
                                      pie = "equal",                   # removed
                                      legend_n = 3,                    # removed
                                      offset_tiplab = rel(1),          # removed
                                      hclust_method = "ward.D",        # removed
                                      group_color = NULL,              # removed
                                      extend = 0.3,                    # removed
                                      hilight = TRUE,                  # removed
                                      geneClusterPanel = "heatMap",    # removed
                                      hexpand = .1,                    # removed
                                      align = "both",                  # removed
                                      cluster.params = list(              
                                          method = "ward.D",     
                                          n = 5,                 
                                          color = NULL,                                        
                                          label_words_n = 4,     
                                          label_format = 30    
                                      ),
                                      hilight.params = list(
                                          hilight = TRUE,                             
                                          align = "both"              
                                      ),                                      
                                      clusterPanel.params = list(
                                          clusterPanel = "heatMap",    
                                          pie = "equal",                
                                          legend_n = 3,
                                          colnames_angle = 0                  
                                      ),
                                      offset.params = list(
                                          bar_tree = rel(1),           
                                          tiplab = rel(1),             
                                          extend = 0.3,                
                                          hexpand = .1                  
                                      ),...) {

    # change parameter name
    ##############################################################
    params_df <- as.data.frame(rbind(
        c("hclust_method", "cluster.params", "method"),
        c("nCluster", "cluster.params", "n"),
        c("group_color", "cluster.params", "color"),
        c("nWords", "cluster.params", "label_words_n"),
        c("label_format_cladelab", "cluster.params", "label_format"),

        c("hilight", "hilight.params", "hilight"),
        c("align", "hilight.params", "align"),

        c("geneClusterPanel", "clusterPanel.params", "clusterPanel"),
        c("pie", "clusterPanel.params", "pie"),
        c("legend_n", "clusterPanel.params", "legend_n"),

        c("offset", "offset.params", "bar_tree"),
        c("offset_tiplab", "offset.params", "tiplab"),
        c("extend", "offset.params", "extend"),
        c("hexpand", "offset.params", "hexpand"))
    )
    colnames(params_df) <- c("original", "listname", "present")
    rownames(params_df) <- params_df$original

    default.cluster.params <- list(
        method = "ward.D",    
        n = 5,                
        color = NULL,                               
        label_words_n = 4,    
        label_format = 30     
    )

    default.hilight.params <- list(
        hilight = TRUE,                       
        align = "both"            
    )
    
    default.clusterPanel.params <- list(
        clusterPanel = "heatMap", 
        pie = "equal",            
        legend_n = 3,
        colnames_angle = 0              
    )

    default.offset.params <- list(
        bar_tree = rel(1),  
        tiplab = rel(1),    
        extend = 0.3,       
        hexpand = .1        
    )
    # use modifyList to change the values of parameter 
    cluster.params <- modifyList(default.cluster.params, cluster.params)   
    hilight.params <- modifyList(default.hilight.params, hilight.params) 
    clusterPanel.params <- modifyList(default.clusterPanel.params, clusterPanel.params) 
    offset.params <- modifyList(default.offset.params, offset.params)  
    params_list <- list(x = x,
        showCategory = showCategory,
        color = color,
        nWords = nWords,                     
        nCluster = nCluster,                   
        cex_category = cex_category, 
        split = split,
        label_format = label_format, 
        label_format_cladelab = label_format_cladelab,     
        label_format_tiplab = label_format_tiplab, 
        fontsize = fontsize, 
        offset = offset,                
        pie = pie,                  
        legend_n = legend_n,                   
        offset_tiplab = offset_tiplab,         
        hclust_method = hclust_method,       
        group_color = group_color,             
        extend = extend,                   
        hilight = hilight,                 
        geneClusterPanel = geneClusterPanel,   
        hexpand = hexpand,                   
        align = align,                 
        cluster.params = cluster.params,
        hilight.params = hilight.params,                                      
        clusterPanel.params = clusterPanel.params,
        offset.params = offset.params
    )     

    # get all parameters value
    args <- as.list(match.call())
    removed_params <- intersect(params_df$original, names(args))
    if (length(removed_params) > 0) {
        for (i in removed_params) {
            params_list[[params_df[i, 2]]][[params_df[i, 3]]] <- get(i)
            warn <- get_param_change_message(i, params_df)
            warning(warn)
        }
    }


    cluster.params <- params_list[["cluster.params"]]  
    hilight.params <- params_list[["hilight.params"]]
    clusterPanel.params <- params_list[["clusterPanel.params"]]
    offset.params <- params_list[["offset.params"]]

    hclust_method <- cluster.params[["method"]]
    nCluster <- cluster.params[["n"]]
    group_color <- cluster.params[["color"]]
    nWords <- cluster.params[["label_words_n"]]
    label_format_cladelab <- cluster.params[["label_format"]]
    hilight <- hilight.params[["hilight"]]
    align <- hilight.params[["align"]]
    geneClusterPanel <- clusterPanel.params[["clusterPanel"]]
    pie <- clusterPanel.params[["pie"]]
    legend_n <- clusterPanel.params[["legend_n"]]
    colnames_angle <- clusterPanel.params[["colnames_angle"]]
    hilight <- hilight.params[["hilight"]]
    align <- hilight.params[["align"]]
    offset <- offset.params[["bar_tree"]]
    offset_tiplab <- offset.params[["tiplab"]]
    extend <- offset.params[["extend"]]
    hexpand <- offset.params[["hexpand"]]

    # geneClusterPanel <- match.arg(geneClusterPanel, c("heatMap", "dotplot", "pie"))                                  
    has_pairsim(x)
    group <- Description <- Cluster <- Count <- NULL
    # to compatible with older versions
    if (!is.null(label_format)) {
        label_format_cladelab <- label_format
    }
    # y <- get_selected_category(showCategory, x, split)  
    y <- fortify(x, showCategory = showCategory,
             includeAll = TRUE, split = split)
    y$Cluster <- sub("\n.*", "", y$Cluster)
    if ("core_enrichment" %in% colnames(y)) { ## for GSEA result
        y$geneID <- y$core_enrichment
    }
    ## Data structure transformation, combining the same ID (Description) genes
    merged_ggData <- merge_compareClusterResult(y)
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    ## Fill the upper triangular matrix completely
    termsim2 <- fill_termsim(x, rownames(ID_Cluster_mat))
    hc <- stats::hclust(stats::as.dist(1- termsim2),
                        method = hclust_method)
    clus <- stats::cutree(hc, nCluster)
    rownames(merged_ggData) <- merged_ggData$Description
    d <- data.frame(label = names(clus),
        count = merged_ggData[names(clus), "Count"])
  
    p <- group_tree(hc = hc, clus = clus, d = d, offset_tiplab = offset_tiplab,
        nWords = nWords, label_format_cladelab = label_format_cladelab, 
        label_format_tiplab = label_format_tiplab, offset = offset, 
        fontsize = fontsize, group_color = group_color, extend = extend, 
        hilight = hilight, cex_category = cex_category, ID_Cluster_mat = ID_Cluster_mat,
        geneClusterPanel = geneClusterPanel, align = align, add_tippoint = FALSE,
        align_tiplab = TRUE, color = color)

     
    p_data <- as.data.frame(p$data)
    p_data <- p_data[which(!is.na(p_data$label)), ]
    rownames(p_data) <- p_data$label
    p_data <- p_data[rownames(ID_Cluster_mat), ]
    # if (geneClusterPanel == "pie") {
    #     xlim <- c(0, xlim * 3 * p_data$x[1])  
    # } else {
    #     xlim <- c(0, xlim * 4.5 * p_data$x[1])
    # }
    
    if (geneClusterPanel == "pie") {
        ID_Cluster_mat$radius <- sqrt(p_data$count / sum(p_data$count) * cex_category)
        ID_Cluster_mat$x <- p_data$x
        ID_Cluster_mat$y <- p_data$y
        ID_Cluster_mat$node <- p_data$node
        p <- p + ggnewscale::new_scale_fill() +
            scatterpie::geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                    cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-4)],color=NA) +
            scatterpie::geom_scatterpie_legend(ID_Cluster_mat$radius,
                x = 0.8, y = 0.1, n = legend_n,
                labeller = function(x) round(sum(p_data$count) * x^2 / cex_category)) +
            # coord_equal(xlim = xlim) +
            labs(fill = "Cluster") + coord_equal()
    }
    
    if (geneClusterPanel == "heatMap") {
        heatMapData <- as.data.frame(x)
        heatMapData$Cluster <- as.character(heatMapData$Cluster)
        heatMapData <- heatMapData[heatMapData$Cluster %in% colnames(ID_Cluster_mat), ]
        heatMapData <- heatMapData[heatMapData$Description %in% rownames(ID_Cluster_mat), ]
        for (i in seq_len(nrow(heatMapData))) {
            ID_Cluster_mat[heatMapData[i, "Description"], heatMapData[i, "Cluster"]] <- heatMapData[i, color]
        }
        
        
        p <- p + ggnewscale::new_scale_fill() # +
            # coord_equal(xlim = xlim)
        p <- ggtree::gheatmap(p, ID_Cluster_mat, colnames_angle = colnames_angle) + 
            # scale_fill_continuous(trans = "log10", name = color) +
            set_enrichplot_color(type = "fill", trans = "log10", name = color)
            
    }

    if (geneClusterPanel == "dotplot") {
        dotdata <- as.data.frame(x)
        pData <- as.data.frame(p$data)
        paths <- pData$label[order(pData$y, decreasing = TRUE)] %>% .[!is.na(.)]
        dotdata <- dotdata[dotdata$Description %in% paths, ]
        dotdata <- dplyr::select(dotdata, Description, dplyr::everything())
        check_installed("ggtreeExtra", "for `treeplot()` with ` clusterPanel = 'dotplot'`.")
	    p <- p + ggnewscale::new_scale_colour() + 
            ggtreeExtra::geom_fruit(data = dotdata, geom = geom_point,
                       mapping = aes_string(x = "Cluster", y = "Description", 
                                     size = "Count", color = color),
                       # pwidth = 0.5, offset = -0.2,
                       pwidth = 0.06*ncol(ID_Cluster_mat),
                       axis.params = list(axis = "x", text.size = 3, line.alpha = 0, text.angle = colnames_angle)) +
            # scale_colour_continuous(trans = "log10", name = color) + 
            set_enrichplot_color(trans = "log10", name = color)
            
    }
    p + ggtree::hexpand(ratio = hexpand)

}




##' Fill the upper triangular matrix completely
##'
##' @param x enrichment result
##' @param keep keep value
##'
##' @return a data.frame
##' @noRd
fill_termsim <- function(x, keep) {
    termsim <- x@termsim[keep, keep]
    termsim[which(is.na(termsim))] <- 0
    termsim2 <- termsim + t(termsim)
    for ( i in seq_len(nrow(termsim2)))
        termsim2[i, i] <- 1
    return(termsim2)
}

##' Add geom_cladelab() to a ggtree object.
##'
##' @param p a ggtree object
##' @param nWords the number of words in the cluster tags
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param offset offset of bar and text from the clade.
##' @importFrom rlang check_installed
##'
##' @return a ggtree object
##' @noRd
add_cladelab <- function(p, nWords, label_format_cladelab, 
                         offset, roots, 
                         fontsize, group_color, cluster_color, 
                         pdata, extend, hilight, align) {
    # align <- getOption("enriplot.treeplot.align", default = "both")
    cluster_label <- sapply(cluster_color, get_wordcloud, ggData = pdata,
                        nWords = nWords)
    label_func_cladelab <- default_labeller(label_format_cladelab)
    if (is.function(label_format_cladelab)) {
        label_func_cladelab <- label_format_cladelab
    }
    cluster_label <- label_func_cladelab(cluster_label)
    n_color <- length(levels(cluster_color)) - length(cluster_color)
    if (is.null(group_color)) {
        check_installed('scales', 'for `add_cladelab()`.')
	color2 <- scales::hue_pal()(length(roots) + n_color)
        if (n_color > 0) color2 <- color2[-seq_len(n_color)]
    } else {
        color2 <- group_color
    }
    df <- data.frame(node = as.numeric(roots),
        labels = cluster_label,
        cluster=cluster_color,
        # color = scales::hue_pal()(length(roots) + n_color)[-seq_len(n_color)]
        color = color2
    )
    
    p <- p + ggnewscale::new_scale_colour() + 
        geom_cladelab(
            data = df,
            mapping = aes_(node =~ node, label =~ labels, color =~ cluster),
            textcolor = "black",
            extend = extend,
            show.legend = FALSE,
            fontsize = fontsize, offset = offset) + 
            scale_color_manual(values = df$color, 
                               guide = 'none')
    if (hilight) {
        p <- p + ggtree::geom_hilight(
            data = df,
            mapping = aes_(node =~ node, fill =~ cluster),
            show.legend = FALSE, 
            align = align) + 
            scale_fill_manual(values = df$color, 
                               guide = 'none')

    }
    
    return(p)
 
}

##' Group the nodes.
##'
##' @return a ggtree object
##' @noRd
group_tree <- function(hc, clus, d, offset_tiplab, nWords, 
                       label_format_cladelab, label_format_tiplab,
                       offset, fontsize, group_color, 
                       extend, hilight, cex_category, 
                       ID_Cluster_mat = NULL, geneClusterPanel = NULL,
                       align, add_tippoint = TRUE, align_tiplab = TRUE, color) {
    group <- count <- NULL
    # cluster data
    dat <- data.frame(name = names(clus), cls=paste0("cluster_", as.numeric(clus)))
    grp <- apply(table(dat), 2, function(x) names(x[x == 1]))  
    p <- ggtree(hc, hang=-1, branch.length = "none", show.legend=FALSE)
    # extract the most recent common ancestor
    noids <- lapply(grp, function(x) unlist(lapply(x, function(i) ggtree::nodeid(p, i))))
    roots <- unlist(lapply(noids, function(x) ggtree::MRCA(p, x)))
    # cluster data
    p <- ggtree::groupOTU(p, grp, "group") + aes_(color =~ group)
    rangeX <- max(p$data$x, na.rm=TRUE) - min(p$data$x, na.rm=TRUE)
    if (inherits(offset_tiplab, "rel")) {
        offset_tiplab <- unclass(offset_tiplab)
        if (geneClusterPanel == "pie" || is.null(geneClusterPanel)) {
        ## 1.5 * max(radius_of_pie)
            offset_tiplab <- offset_tiplab * .5 * max(sqrt(d$count / sum(d$count) * cex_category))
        }  else if (geneClusterPanel == "heatMap") {
            ## Close to the width of the tree
            offset_tiplab <- offset_tiplab * 0.2 * ncol(ID_Cluster_mat) * rangeX
        } else if (geneClusterPanel == "dotplot") {
            ## Close to the width of the tree
            offset_tiplab <- offset_tiplab * 0.09 * ncol(ID_Cluster_mat) * rangeX
        }
    }

    if (inherits(offset, "rel")) {
        offset <- unclass(offset)
        offset <- offset * rangeX * 1.2 + offset_tiplab 
    }
    # max_nchar <- max(nchar(p$data$label), na.rm = TRUE)
       
    pdata <- data.frame(name = p$data$label, color2 = p$data$group)
    pdata <- pdata[!is.na(pdata$name), ]
    cluster_color <- unique(pdata$color2)
    n_color <- length(levels(cluster_color)) - length(cluster_color)
    if (!is.null(group_color)) {
        color2 <- c(rep("black", n_color), group_color)
        p <- p + scale_color_manual(values = color2, guide = 'none')
    }
    p <- p %<+% d


    if (!is.null(label_format_tiplab)) {
        label_func_tiplab <- default_labeller(label_format_tiplab)
        if (is.function(label_format_tiplab)) {
            label_func_tiplab <- label_format_tiplab
        }
        isTip <- p$data$isTip
        p$data$label[isTip] <-  label_func_tiplab(p$data$label[isTip])
    }

    p <- add_cladelab(p = p, nWords = nWords, 
        label_format_cladelab = label_format_cladelab,
        offset = offset, roots = roots, fontsize = fontsize, 
        group_color = group_color, cluster_color = cluster_color, 
        pdata = pdata, extend = extend, hilight = hilight, align = align)
    if (add_tippoint) {
        p <- p + ggnewscale::new_scale_colour() +
            geom_tippoint(aes(color = color, size = count)) + 
            # scale_colour_continuous(name = color)+
            set_enrichplot_color(name = color)
    }
    ## add tiplab 
    p <- p + geom_tiplab(offset = offset_tiplab, hjust = 0,
                show.legend = FALSE, align = align_tiplab, linesize = 0)     
    return(p)
}
