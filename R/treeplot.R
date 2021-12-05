##' @rdname treeplot
##' @exportMethod treeplot
setMethod("treeplot", signature(x = "enrichResult"),
    function(x, showCategory = 30, color = "p.adjust", label_format = 30, ...) {
        treeplot.enrichResult(x, showCategory = showCategory,
            color = color, label_format = label_format, ...)
    })

##' @rdname treeplot
##' @exportMethod treeplot
setMethod("treeplot", signature(x = "gseaResult"),
    function(x, showCategory = 30, color = "p.adjust", label_format = 30, ...) {
        treeplot.enrichResult(x, showCategory = showCategory,
            color = color, label_format = label_format, ...)
    })

##' @rdname treeplot
##' @exportMethod treeplot
setMethod("treeplot", signature(x = "compareClusterResult"),
    function(x, showCategory = 5, color = "p.adjust", label_format = 30, ...) {
        treeplot.compareClusterResult(x, showCategory = showCategory,
            color = color, label_format = label_format, ...)
    })



##' @rdname treeplot
##' @param nWords The number of words in the cluster tags.
##' @param nCluster The number of clusters, the default value is 5.
##' @param cex_category Number indicating the amount by which plotting category.
##' nodes should be scaled relative to the default.
## ' @param xlim Limits for the x axes, the default value is 1. If the picture is not 
##' displayed completely, the user can increase this value.
##' @param offset_tiplab tiplab offset, the bigger the number, 
##' the farther the distance between the node and the branch. 
##' The default is 1, when geneClusterPanel = "pie", meaning 1 *  max_radius_of_the_pies; 
##' when geneClusterPanel = "heatMap", meaning 1 * 0.16 * column_number_of_heatMap * x_range_of_tree;
##' when geneClusterPanel = "dotplot", meaning 1 * 0.09 * column_number_of_dotplot * x_range_of_tree.
##' @param offset numeric, distance bar and tree, offset of bar and text from the clade, default is 1,
##' meaning 1 * 1.2 * x_range_of_tree plus distance_between_tree_and_tiplab
##' (1 * (1.2 * x_range_of_tree + distance_between_tree_and_tiplab)).
##' @param fontsize The size of text, default is 4.
##' @param hclust_method Method of hclust. This should be (an unambiguous abbreviation of) one of "ward.D", 
##' "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
##' @param group_color A vector of group colors, the length of the vector should be the same as nCluster.
##' @param extend Numeric, extend the length of bar, default is 0.3.
##' @param hilight Logical value, if TRUE(default), add ggtree::geom_hilight() layer.
##' @param hexpand expand x limits by amount of xrange * hexpand.
##' @param align control the align direction of the edge of high light rectangular.
##' Options is 'none', 'left', 'right', 'both (default)'.
##' @importFrom ggtree `%<+%`
##' @importFrom ggtree ggtree
##' @importFrom ggtree geom_tiplab
##' @importFrom ggtree geom_tippoint
##' @importFrom ggtree groupClade
##' @importFrom ggtree geom_cladelab
##' @importFrom ggplot2 coord_cartesian
##' @importFrom ggplot2 scale_colour_continuous
##'
treeplot.enrichResult <- function(x, showCategory = 30,
                                  color = "p.adjust",
                                  nWords = 4, nCluster = 5,
                                  cex_category = 1,
                                  label_format = 30, # xlim = 1,
                                  fontsize = 4, offset = 1,
                                  offset_tiplab = 1, 
                                  hclust_method = "ward.D", 
                                  group_color = NULL, 
                                  extend = 0.3, hilight = TRUE, 
                                  hexpand = .1, align = "both", ...) {
    group <- p.adjust <- count<- NULL

    if (class(x) == "gseaResult")
        x@result$Count <- x$core_enrichment %>%
            strsplit(split = "/")  %>%
            vapply(length, FUN.VALUE = 1)   
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
        nWords = nWords, label_format = label_format, offset = offset, 
        fontsize = fontsize, group_color = group_color, extend = extend, 
        hilight = hilight, cex_category = cex_category, align = align)
    # xlim <-  c(0, xlim * 3 * max(p$data$x))
    # p + coord_cartesian(xlim = xlim) +
      p + ggnewscale::new_scale_colour() +
        geom_tippoint(aes(color = color, size = count)) +
        scale_colour_continuous(low="red", high="blue", name = color, 
            guide = guide_colorbar(reverse = TRUE)) +
        scale_size_continuous(name = "number of genes",
                              range = c(3, 8) * cex_category) + 
        ggtree::hexpand(ratio = hexpand) + 
        guides(size  = guide_legend(order = 1), 
               color = guide_colorbar(order = 2))
}


##' @rdname treeplot
##' @param split Separate result by 'category' variable.
##' @param legend_n Number of circle in legend, the default value is 3.
##' @param geneClusterPanel one of "heatMap"(default), "dotplot", "pie".
##' @param pie Used only when geneClusterPanel = "pie", 
##' proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'.
##' @importFrom ggtree nodepie
##' @importFrom ggtree geom_inset
##' @importFrom ggplot2 scale_fill_manual
treeplot.compareClusterResult <-  function(x, showCategory = 5,
                                      color = "p.adjust",
                                      nWords = 4, nCluster = 5,
                                      cex_category = 1, split = NULL,
                                      label_format = 30, # xlim = 1,
                                      fontsize = 4, offset = 1, pie = "equal",
                                      legend_n = 3, offset_tiplab = 1, 
                                      hclust_method = "ward.D", group_color = NULL, 
                                      extend = 0.3, hilight = TRUE, 
                                      geneClusterPanel = "heatMap", 
                                      hexpand = .1, align = "both", ...) {
    geneClusterPanel <- match.arg(geneClusterPanel, c("heatMap", "dotplot", "pie"))                                  
    has_pairsim(x)
    group <- Description <- Cluster <- Count <- NULL
    # y <- get_selected_category(showCategory, x, split)  
    y <- fortify(x, showCategory = showCategory,
             includeAll = TRUE, split = split)
    y$Cluster <- sub("\n.*", "", y$Cluster)
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
        nWords = nWords, label_format = label_format, offset = offset, 
        fontsize = fontsize, group_color = group_color, extend = extend, 
        hilight = hilight, cex_category = cex_category, ID_Cluster_mat = ID_Cluster_mat,
        geneClusterPanel = geneClusterPanel, align = align)
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
            labs(fill = "Cluster")
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
        p <- ggtree::gheatmap(p, ID_Cluster_mat) + 
            scale_fill_continuous(low="red", high="blue", 
                                  guide = guide_colorbar(reverse=TRUE),
                                  trans = "log10",
                                  name = color) 
    }

    if (geneClusterPanel == "dotplot") {
        dotdata <- as.data.frame(x)
        pData <- as.data.frame(p$data)
        paths <- pData$label[order(pData$y, decreasing = TRUE)] %>% .[!is.na(.)]
        dotdata <- dotdata[dotdata$Description %in% paths, ]
        dotdata <- dplyr::select(dotdata, Description, dplyr::everything())
        p <- p + ggnewscale::new_scale_colour() + 
            ggtreeExtra::geom_fruit(data = dotdata, geom = geom_point,
                       mapping = aes_string(x = "Cluster", y = "Description", 
                                     size = "Count", color = color),
                       pwidth = 0.5,
                       axis.params = list(axis = "x", text.size = 3, line.alpha = 0)) +
            scale_colour_continuous(low="red", high="blue", 
                                  guide=guide_colorbar(reverse=TRUE),
                                  trans = "log10",
                                  name = color) # + 
            # coord_equal(xlim = xlim) 
            
    }
    p + ggtree::hexpand(ratio = hexpand) + coord_equal()
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
##'
##' @return a ggtree object
##' @noRd
add_cladelab <- function(p, nWords, label_format, offset, roots, 
                         fontsize, group_color, cluster_color, 
                         pdata, extend, hilight, align) {
    # align <- getOption("enriplot.treeplot.align", default = "both")
    cluster_label <- sapply(cluster_color, get_wordcloud, ggData = pdata,
                        nWords = nWords)
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    cluster_label <- label_func(cluster_label)
    #names(cluster_label) <- cluster_color
    n_color <- length(levels(cluster_color)) - length(cluster_color)
    if (is.null(group_color)) {
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
                       label_format, offset, fontsize, group_color, 
                       extend, hilight, cex_category, 
                       ID_Cluster_mat = NULL, geneClusterPanel = NULL,
                       align) {
    group <- NULL
    # cluster data
    dat <- data.frame(name = names(clus), cls=paste0("cluster_", as.numeric(clus)))
    grp <- apply(table(dat), 2, function(x) names(x[x == 1]))  
    p <- ggtree(hc, branch.length = "none", show.legend=FALSE)
    # extract the most recent common ancestor
    noids <- lapply(grp, function(x) unlist(lapply(x, function(i) ggtree::nodeid(p, i))))
    roots <- unlist(lapply(noids, function(x) ggtree::MRCA(p, x)))
    # cluster data
    p <- ggtree::groupOTU(p, grp, "group") + aes_(color =~ group)


    if (geneClusterPanel == "pie" || is.null(geneClusterPanel)) {
        ## 1.5 * max(radius_of_pie)
        offset_tiplab <- offset_tiplab * 1.5 * max(sqrt(d$count / sum(d$count) * cex_category))
    }  else if (geneClusterPanel == "heatMap") {
        ## Close to the width of the tree
        offset_tiplab <- offset_tiplab * 0.16 * ncol(ID_Cluster_mat) * max(p$data$x)
    } else if (geneClusterPanel == "dotplot") {
        ## Close to the width of the tree
        offset_tiplab <- offset_tiplab * 0.09 * ncol(ID_Cluster_mat) * max(p$data$x)
    }
    # max_nchar <- max(nchar(p$data$label), na.rm = TRUE)
    offset <- offset * (max(p$data$x) * 1.2 + offset_tiplab)    
    pdata <- data.frame(name = p$data$label, color2 = p$data$group)
    pdata <- pdata[!is.na(pdata$name), ]
    cluster_color <- unique(pdata$color2)
    n_color <- length(levels(cluster_color)) - length(cluster_color)
    if (!is.null(group_color)) {
        color2 <- c(rep("black", n_color), group_color)
        p <- p + scale_color_manual(values = color2, guide = 'none')
    }
    
    p <- p %<+% d +
        geom_tiplab(offset = offset_tiplab, hjust = 0,
                    show.legend = FALSE, align = TRUE, linesize = 0)
    
    p <- add_cladelab(p = p, nWords = nWords, label_format = label_format, 
        offset = offset, roots = roots, fontsize = fontsize, 
        group_color = group_color, cluster_color = cluster_color, 
        pdata = pdata, extend = extend, hilight = hilight, align = align) 
    return(p)
}
