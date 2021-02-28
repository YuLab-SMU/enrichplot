##' @rdname mdsplot
##' @exportMethod mdsplot
setMethod("mdsplot", signature(x = "enrichResult"),
    function(x, showCategory = 30, color = "p.adjust", label_format = 30, ...) {
        mdsplot.enrichResult(x, showCategory = showCategory,
            color = color, label_format = label_format, ...)
    })

##' @rdname mdsplot
##' @exportMethod mdsplot
setMethod("mdsplot", signature(x = "gseaResult"),
    function(x, showCategory = 30, color = "p.adjust", label_format = 30, ...) {
        mdsplot.enrichResult(x, showCategory = showCategory,
            color = color, label_format = label_format, ...)
    })

##' @rdname mdsplot
##' @exportMethod mdsplot
setMethod("mdsplot", signature(x = "compareClusterResult"),
    function(x, showCategory = 30, label_format = 30, ...) {
        mdsplot.compareClusterResult(x, showCategory = showCategory,
            label_format = label_format, ...)
    })





##' @rdname mdsplot
##' @param x Enrichment result.
##' @param showCategory A number or a vector of terms. If it is a number, 
##' the first n terms will be displayed. If it is a vector of terms, 
##' the selected terms will be displayed.
##' @param color variable that used to color enriched terms, e.g. pvalue,
##' p.adjust or qvalue
##' @param nWords Numeric, the number of words in the cluster tags, the default value is 4.
##' @param nCluster Numeric, the number of clusters,
##' the default value is square root of the number of nodes.
##' @param split Separate result by 'category' variable
##' @param cex_label_group Numeric, scale of group labels size, the default value is 1.
##' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
##' The default is FALSE
##' @param cex_category Numeric, indicating the size by which plotting category
##' nodes should be scaled relative to the default.
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param repel whether to correct the position of the label. Defaults to FALSE.
##' @param shadowtext a logical value, whether to use shadow font. Defaults to TRUE.
##' @param engine The function used for MDS, 
##' one of "cmdscale" (the default), "sammon", "monoMDS" and "isoMDS".
##' @param group_category Logical, if TRUE (the default), group the terms.
##' @param cex_label_category Scale of category node label size.
##' @param ... Additional parameters used to set the position of the group label.
##' When the parameter repel is set to TRUE, additional parameters will take effect.
##'
##' additional parameters can refer the following parameters.
##'     \itemize{
##'        \item \code{force} Force of repulsion between overlapping text labels. Defaults to 1.
##'        \item \code{nudge_x, nudge_y} Horizontal and vertical adjustments to nudge
##'         the starting position of each text label.
##'        \item \code{direction} "both", "x", or "y" – direction in which to adjust position of labels.
##'     }
##'
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_size_continuous
##' @importFrom stats kmeans
##' @export
mdsplot.enrichResult <- function(x, showCategory = 30,
                    color = "p.adjust", label_format = 30,
                    nWords = 4, nCluster = NULL,
                    split = NULL,
                    cex_label_group = 1,
                    group_legend = FALSE, cex_category = 1,
                    repel = FALSE,
                    shadowtext = TRUE, 
                    engine = "cmdscale", 
                    group_category = TRUE, 
                    cex_label_category = 1, ...) {
    has_pairsim(x)
    label_group <- 3
    label_size_category <- 3
    y <- NULL
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

    termsim2 <- fill_termsim(x, keep)
    data <- get_cluster(termsim2, nCluster = nCluster, engine = engine)
    terms <- rownames(termsim2)
    object <- as.data.frame(x)
    rownames(object) <- object$Description
    data <- cbind(data, object[terms, ])

    p <- ggplot(data, aes(x=x, y=y))+
        geom_point(aes_string(size = "Count", fill = color), shape = 21) +
        scale_size_continuous(name = "number of genes",
                              range = c(3, 8) * cex_category) +
        scale_fill_continuous(low = "red", high = "blue", name = color,
                              guide = guide_colorbar(reverse = TRUE)) +
        theme_classic()
    p$data$name <- rownames(p$data)
    if (!group_category) 
        p <-  add_node_label(p = p, data = NULL, label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext)

    add_ellipse(group_legend = group_legend, p = p, data = data, nWords = nWords, 
        label_format = label_format, repel = repel, shadowtext = shadowtext, 
        label_group = label_group, cex_label_group = cex_label_group, 
        group_category = group_category)
}


##' @rdname mdsplot
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 coord_equal
##' @importFrom scatterpie geom_scatterpie
##' @importFrom scatterpie geom_scatterpie_legend
##' @param pie Proportion of clusters in the pie chart, one of
##' 'equal' (default) and 'Count'.
##' @param legend_n Number of circle in legend, the default value is 5.
mdsplot.compareClusterResult <- function(x, showCategory = 30,
                                         nWords = 4, nCluster = NULL, 
                                         split = NULL, 
                                         cex_label_group = 1, 
                                         pie = "equal", legend_n = 3,
                                         cex_category = 1, 
                                         group_legend = FALSE,
                                         label_format = 30, 
                                         repel = FALSE, shadowtext = TRUE,
                                         engine = "cmdscale", 
                                         group_category = TRUE, 
                                         cex_label_category = 1, ...) {
    has_pairsim(x)
    label_group <- 3
    label_size_category <- 3
    y <- get_selected_category(showCategory, x, split)
    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- merge_compareClusterResult(y)
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    ## Fill the upper triangular matrix completely

    termsim2 <- fill_termsim(x = x, keep = rownames(ID_Cluster_mat))
    data <- get_cluster(termsim2  = termsim2, nCluster = nCluster, engine = engine)
    data <- cbind(data, y_union)
    ID_Cluster_mat$radius <- sqrt(data$Count / sum(data$Count) * cex_category * 0.0125)
    ID_Cluster_mat$x <- data$x
    ID_Cluster_mat$y <- data$y
    x_loc1 <- min(ID_Cluster_mat$x)
    y_loc1 <- min(ID_Cluster_mat$y)
    p <- ggplot(data, aes(x=x, y=y)) +
        scatterpie::geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
            cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-4)],color=NA) +
        coord_equal()+
        theme_classic() +
        scatterpie::geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
            n = legend_n,
            labeller = function(x) round(sum(data$Count) * x^2 / (cex_category * 0.0125))) +
        theme(legend.title = element_text(size = 10),
              legend.text  = element_text(size = 10)) +
        theme(panel.background = element_blank())
    p$data$name <- rownames(p$data)
    if (!group_category) 
        p <-  add_node_label(p = p, data = NULL, label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext)

    add_ellipse(group_legend = group_legend, p = p, data = data, nWords = nWords, 
        label_format = label_format, repel = repel, shadowtext = shadowtext, 
        label_group = label_group, cex_label_group = cex_label_group, 
        group_category = group_category)
}





#' Group terms
#'
#' @param termsim2 Similarity matrix of terms.
#' @param nCluster Numeric, the number of cluster.
#'
#' @noRd
get_cluster <- function(termsim2 = termsim2, nCluster, engine) {
    ## If the similarity between the two terms is 1, an error will be reported, so fine-tuning.
    termsim2[which(termsim2 == 1)] <- 0.99999
    for (i in seq_len(nrow(termsim2))) termsim2[i, i] <- 1
    
    dune.dist <- stats::as.dist(1- termsim2)
    # res <- ape::pcoa(dune.dist, correction="none")
    if (engine == "cmdscale") {
        res <- stats::cmdscale(dune.dist, eig = T)   
    }

    if (engine == "sammon") {
        y <- stats::cmdscale(dune.dist)
        ## If the matrix y has duplicate rows it will report an error, so perturb slightly  
        dup <- which(duplicated(y) == TRUE)  
        y[dup, 1] <- y[dup, 1] + 10^-7 * seq_len(length(dup))  
        res <- MASS::sammon(d = dune.dist, y = y)
    }  
    
    if (engine == "monoMDS" || engine == "isoMDS") {
        res <- vegan::metaMDS(dune.dist, engine = engine)
    }
    data <- as.data.frame(res$points[, 1:2])
    colnames(data) <- c('x', 'y')

    if (is.null(nCluster)){
        data$color <- kmeans(data, floor(sqrt(nrow(data))))$cluster
    } else {
        if (nCluster > nrow(data)) nCluster <- nrow(data)
        data$color <- kmeans(data, nCluster)$cluster
    }
    data$color <- as.character(data$color)
    if (engine == "cmdscale") {
        data$pocas <- 0
        data$pocas[seq_len(length(res$eig))] <- as.numeric(res$eig)
    } else {
        data$stress <- res$stress
    }
    return(data)
}

##' Add circle and group label to a plot.
##'
##' @param group_legend A logical, if TRUE, add the group legend.
##' @param p A ggplot2 object
##' @param data a data.frame containing group information
##' @param nWords Numeric, the number of words in the cluster tags.
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param repel Whether to correct the position of the label. Defaults to FALSE.
##' @param shadowtext A logical value, whether to use shadow font.
##' @param label_group The basic size of gourp label.
##' @param cex_label_group Numeric, scale of group labels size.
##'
##' @noRd
add_ellipse <- function(group_legend, p, data, nWords,
                        label_format, repel, shadowtext, label_group,
                        cex_label_group, group_category) {
    pocas <- data$pocas
    show_legend <- c(group_legend, FALSE)
    names(show_legend) <- c("fill", "color")
    if ("pocas" %in% colnames(data)) {
        xlab = paste("pcoa 1 (", format(100 * pocas[1] / sum(pocas), digits=4), "%)", sep = "")
        ylab = paste("pcoa 2 (", format(100 * pocas[2] / sum(pocas), digits=4), "%)", sep = "")
        title = "PCoA"
    } else {
        xlab = "NMDS1"
        ylab = "NMDS2"
        title = paste("NMDS (stress = ", format(data$stress[1], digits=4), ")", sep = "")
    }
    p <- p +  labs(x = xlab, y = ylab, title = title) 
    if (group_category) {
        p <- p + ggnewscale::new_scale_fill() +
            ggforce::geom_mark_ellipse(aes_(x = ~x, y = ~y, color = ~color,
                         fill =~ color), show.legend = show_legend) +
            theme(legend.title = element_text(size = 10),
                       legend.text  = element_text(size = 10)) +
            theme(panel.background = element_blank())
        goid <- data$ID
        data$name <- data$Description
        cluster_color <- unique(data$color)
        clusters <- lapply(cluster_color, function(i){goid[which(data$color == i)]})
        cluster_label <- sapply(cluster_color, get_wordcloud, pdata2 = data,
                                nWords = nWords)
        names(cluster_label) <- cluster_color
        data$color <- cluster_label[data$color]
        label_location <- get_label_location(data, label_format)
        p <- add_group_label(repel = repel, shadowtext = shadowtext, p = p,
            label_location = label_location, label_group = label_group,
            cex_label_group = cex_label_group)
    }
    return(p) 
}


