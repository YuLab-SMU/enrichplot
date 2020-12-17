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
##' @param nWords the number of words in the cluster tags
##' @param nCluster the number of clusters
##' @param cex_category number indicating the amount by which plotting category nodes should be scaled relative to the default.
##' @param xlim Limits for the x axes
##' @param offset distance bar and tree, offset of bar and text from the clade, default is 0.9.
##' @param fontsize the size of text, default is 4.
##' @param offset_tiplab tiplab offset
##' @param hclust_method method of hclust. This should be (an unambiguous abbreviation of) one of "ward.D", 
##' "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
##' @param group_color a vector of group colors.
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
                                  label_format = 30, xlim = NULL,
                                  fontsize = 4, offset = NULL,
                                  offset_tiplab = 0.1, 
                                  hclust_method = "ward.D", group_color = NULL, ...) {
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

    ## Group the nodes and find the root node of each group of nodes.
    p <- group_tree(hc, clus, d, offset_tiplab, nWords, 
        label_format, offset, fontsize, group_color)
    if(is.null(xlim)) xlim <- c(0, 3 * p$data$x[1])
    p + coord_cartesian(xlim = xlim) +
        ggnewscale::new_scale_colour() +
        geom_tippoint(aes(color = color, size = count)) +
        scale_colour_continuous(name = color, guide = guide_colorbar(reverse = TRUE)) +
        scale_size_continuous(name = "number of genes",
                              range = c(3, 8) * cex_category)
}


##' @rdname treeplot
##' @param pie proportion of clusters in the pie chart, one of
##' 'equal' (default) or 'Count'
##' @param split separate result by 'category' variable
##' @param legend_n number of circle in legend
##' @importFrom ggtree nodepie
##' @importFrom ggtree geom_inset
##' @importFrom ggplot2 scale_fill_manual
treeplot.compareClusterResult <-  function(x, showCategory = 5,
                                      color = "p.adjust",
                                      nWords = 4, nCluster = 5,
                                      cex_category = 1, split = NULL,
                                      label_format = 30, xlim = NULL,
                                      fontsize = 4, offset = NULL, pie = "equal",
                                      legend_n = 3, offset_tiplab = 0.5, 
                                      hclust_method = "ward.D", group_color = NULL, ...) {
    group <- NULL
    if (is.numeric(showCategory)) {
        y <- fortify(x, showCategory = showCategory,
                                      includeAll = TRUE, split = split)
        y_union <- merge_compareClusterResult(y)
    } else {
        y <- fortify(x, showCategory=NULL,
                                      includeAll = TRUE, split = split)
        n <- update_n(y_union, showCategory)
        y_union <- merge_compareClusterResult(y)
        y_union <- y_union[match(n, y_union$Description),]
    }

    y$Cluster <- sub("\n.*", "", y$Cluster)
    y <- y[y$ID %in% y_union$ID, ]
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    ## Fill the upper triangular matrix completely
    termsim2 <- fill_termsim(x, rownames(ID_Cluster_mat))
    hc <- stats::hclust(stats::as.dist(1- termsim2),
                        method = hclust_method)
    clus <- stats::cutree(hc, nCluster)
    rownames(y_union) <- y_union$Description
    d <- data.frame(label = names(clus),
        count = y_union[names(clus), "Count"])
        
    p <- group_tree(hc, clus, d, offset_tiplab, nWords, 
        label_format, offset, fontsize, group_color)
    p_data <- as.data.frame(p$data)
    p_data <- p_data[which(!is.na(p_data$label)), ]
    rownames(p_data) <- p_data$label
    p_data <- p_data[rownames(ID_Cluster_mat), ]

    ID_Cluster_mat$radius <- sqrt(p_data$count / sum(p_data$count) * cex_category)
    ID_Cluster_mat$x <- p_data$x
    ID_Cluster_mat$y <- p_data$y
    ID_Cluster_mat$node <- p_data$node
    if(is.null(xlim)) xlim <- c(0, 5 * p_data$x[1])

    p + ggnewscale::new_scale_colour() +
        scatterpie::geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-4)],color=NA) +
        scatterpie::geom_scatterpie_legend(cex_category * ID_Cluster_mat$radius,
        x = 0.8, y = 0.1, n = legend_n,
        labeller = function(x) round(sum(p_data$count) * x^2 / cex_category)) +
        coord_equal(xlim = xlim) +
        labs(fill = "Cluster")
}



##' Get all parent nodes of a vectory of nodes
##'
##' @param tree_data a phylo object
##' @param nodes a vector of nodes
##'
##' @return a list of parent nodes
##'
##' @noRd
get_parents <- function(tree_data, nodes) {
    lapply(nodes, function(x) tidytree::ancestor(.data = tree_data, .node = x))
}


##' Title Fill the upper triangular matrix completely
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

add_cladelab <- function(p, nWords, label_format, offset, roots, 
                         fontsize, group_color, cluster_color, pdata) {
    #node <- label <- color <- NULL
    # pdata <- data.frame(name = p$data$label, color = p$data$group)
    # pdata <- pdata[!is.na(pdata$name), ]
    # cluster_color <- unique(pdata$color)
    cluster_label <- sapply(cluster_color, wordcloud_i, pdata2 = pdata,
                        nWords = nWords)
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    cluster_label <- label_func(cluster_label)
    #names(cluster_label) <- cluster_color
    if(is.null(offset)) offset <- p$data$x[1]
    n_color <- length(levels(cluster_color)) - length(cluster_color)
    if (is.null(group_color)) {
        color2 <- scales::hue_pal()(length(roots) + n_color)[-seq_len(n_color)]
    } else {
        color2 <- group_color
    }
    df <- data.frame(node = as.numeric(roots),
        label = cluster_label,
        cluster=cluster_color,
        # color = scales::hue_pal()(length(roots) + n_color)[-seq_len(n_color)]
        color = color2
    )
    
    p <- p + ggnewscale::new_scale_colour() + 
        geom_cladelab(
            data = df,
            mapping = aes_(node =~ node, label =~ label, color =~ cluster),
            textcolor = "black",
            extend = 0.3,
            show.legend = F,
            fontsize = fontsize, offset = offset) + 
            scale_color_manual(values = df$color, # limits = df$color, 
                               guide = FALSE)
    return(p)
 
}

group_tree <- function(hc, clus, d, offset_tiplab, nWords, 
                       label_format, offset, fontsize, group_color) {
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
    pdata <- data.frame(name = p$data$label, color = p$data$group)
    pdata <- pdata[!is.na(pdata$name), ]
    cluster_color <- unique(pdata$color)
    n_color <- length(levels(cluster_color)) - length(cluster_color)
    if (!is.null(group_color)) {
        color2 <- c(rep("black", n_color), group_color)
        p <- p + scale_color_manual(values = color2, guide = FALSE)
    }
    
    p <- p %<+% d +
        geom_tiplab(offset = offset_tiplab, hjust = 0, show.legend = FALSE, align=TRUE)

    p <- add_cladelab(p, nWords, label_format, offset, roots, 
        fontsize, group_color, cluster_color, pdata) 
    return(p)
}
