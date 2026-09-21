#' Tree plot for enrichment results
#'
#' Creates hierarchical tree visualization of enriched terms based on similarity
#'
#' @rdname treeplot
#' @exportMethod treeplot
setMethod("treeplot", signature(x = "enrichResult"), function(x, ...) {
    treeplot_internal(x, size_var = "Count", ...)
})


#' @rdname treeplot
#' @exportMethod treeplot
setMethod("treeplot", signature(x = "gseaResult"), function(x, ...) {
    treeplot_internal(x, size_var = "setSize", ...)
})

#' @rdname treeplot
#' @exportMethod treeplot
setMethod("treeplot", signature(x = "compareClusterResult"), function(x, ...) {
    treeplot_compareCluster(x, ...)
})


#' @rdname treeplot
#' @param showCategory number of enriched terms to display
#' @param color variable to color nodes, e.g. 'p.adjust', 'pvalue', or 'qvalue'
#' @param size_var variable for node size, e.g. 'Count' (for enrichResult) or 'setSize' (for gseaResult)
#' @param nCluster number of clusters for tree cutting
#' @param cluster_method hierarchical clustering method
#' @param label_format wrap length for labels or custom formatting function
#' @param fontsize_tiplab font size for tip labels
#' @param fontsize_cladelab font size for clade labels
#' @param group_color vector of colors for groups
#' @param extend extend length for clade labels
#' @param hilight whether to highlight clades
#' @param align alignment for highlight rectangles
#' @param hexpand expand x limits by amount of xrange * hexpand
#' @param tiplab_offset offset for tip labels
#' @param cladelab_offset offset for clade labels
#' @return ggplot2 object representing the tree plot
#' @importFrom ggtree ggtree geom_tiplab geom_tippoint groupClade geom_cladelab geom_hilight
#' @importFrom ggplot2 scale_size_continuous guides guide_legend guide_colorbar
#' @importFrom stats hclust cutree as.dist
treeplot_internal <- function(
    x,
    showCategory = 30,
    color = "p.adjust",
    size_var = c("Count", "setSize"),
    nCluster = 5,
    cluster_method = "ward.D",
    label_format = 30,
    fontsize_tiplab = 4,
    fontsize_cladelab = 4,
    group_color = NULL,
    extend = 0.3,
    hilight = TRUE,
    align = "both",
    hexpand = 0.1,
    tiplab_offset = 0.2,
    cladelab_offset = 1
) {
    # Input validation
    if (!inherits(x, c("enrichResult", "gseaResult", "compareClusterResult"))) {
        stop(
            "x must be an enrichResult, gseaResult, or compareClusterResult object"
        )
    }

    # Get selected categories
    n <- update_n(x, showCategory)
    if (is.numeric(n)) {
        keep <- seq_len(n)
    } else {
        keep <- match(n, rownames(x@termsim))
    }

    if (length(keep) == 0) {
        stop("no enriched term found...")
    }

    # Prepare similarity matrix
    termsim2 <- fill_termsim(x, keep)

    # Single-pathway boundary: hierarchical clustering needs at least two objects.
    if (length(keep) == 1) {
        d <- data.frame(
            label = rownames(termsim2),
            stringsAsFactors = FALSE
        )
        return(
            ggplot(d, aes(x = 0, y = 0, label = .data$label)) +
                geom_text() +
                theme_void()
        )
    }

    # Hierarchical clustering
    hc <- hclust(as.dist(1 - termsim2), method = cluster_method)
    nCluster <- min(nCluster, max(1, length(keep)))
    clus <- cutree(hc, nCluster)

    # Prepare data for plotting
    size_var <- intersect(size_var, colnames(x[]))[1]
    if (is.na(size_var)) {
        stop("size_var not found in enrichment result")
    }
    
    # Extract columns and ensure they exist
    d <- x[keep, c(color, size_var)]
    
    # Handle case where columns are collapsed (e.g. tibble with duplicate columns)
    if (ncol(d) == 1 && color == size_var) {
        d <- data.frame(d, d)
        names(d) <- c(color, size_var)
    }

    # Determine safe size column name
    size_col <- "size"
    if (color == "size") {
        size_col <- "size_value"
    }

    # Rename size column
    names(d)[2] <- size_col
    
    # Add label column from cluster names
    d$label <- names(clus)
    
    # Select columns safely
    d <- d[, c("label", color, size_col)]

    # Create tree plot
    p <- create_tree_plot(
        hc = hc,
        clus = clus,
        data = d,
        label_format = label_format,
        fontsize_tiplab = fontsize_tiplab,
        fontsize_cladelab = fontsize_cladelab,
        group_color = group_color,
        extend = extend,
        hilight = hilight,
        align = align,
        color_var = color,
        size_var = size_col,
        tiplab_offset = tiplab_offset,
        cladelab_offset = cladelab_offset
    )

    # Add styling
    p <- p +
        scale_size_continuous(
            name = size_var,
            range = c(3, 8)
        ) +
        ggtree::hexpand(ratio = hexpand) +
        guides(
            size = guide_legend(order = 1),
            color = guide_colorbar(order = 2)
        )

    return(p)
}

#' Tree plot for compareClusterResult objects
#'
#' @param x compareClusterResult object
#' @param showCategory number of enriched terms to display
#' @param color variable to color nodes
#' @param nCluster number of clusters
#' @param cluster_method hierarchical clustering method
#' @param label_format label formatting
#' @param fontsize_tiplab tip label font size
#' @param fontsize_cladelab clade label font size
#' @param group_color group colors
#' @param extend extend length
#' @param hilight whether to highlight clades
#' @param align highlight alignment
#' @param hexpand expand x limits
#' @param tiplab_offset tip label offset
#' @param cladelab_offset clade label offset
#' @param pie proportion method for pie charts ("equal" or "Count")
#' @param cluster_panel panel type for clusters ("pie", "heatMap", or "dotplot")
#' @param legend_n number of legend items for pie charts
#' @param colnames_angle angle for column names in heatmaps
#' @importFrom ggtree gheatmap
#' @importFrom scatterpie geom_scatterpie geom_scatterpie_legend
#' @importFrom ggnewscale new_scale_fill new_scale_colour
#' @noRd
treeplot_compareCluster <- function(
    x,
    showCategory = 30,
    color = "p.adjust",
    nCluster = 5,
    cluster_method = "ward.D",
    label_format = 30,
    fontsize_tiplab = 4,
    fontsize_cladelab = 4,
    group_color = NULL,
    extend = 0.3,
    hilight = TRUE,
    align = "both",
    hexpand = 0.1,
    tiplab_offset = 0.2,
    cladelab_offset = 1,
    pie = "equal",
    cluster_panel = "pie",
    legend_n = 3,
    colnames_angle = 0
) {
    # Prepare data for compareClusterResult
    y <- fortify(
        x,
        showCategory = showCategory,
        includeAll = TRUE,
        split = NULL
    )
    y$Cluster <- sub("\n.*", "", y$Cluster)

    if ("core_enrichment" %in% colnames(y)) {
        y$geneID <- y$core_enrichment
    }

    # Prepare cluster matrix
    ID_Cluster_mat <- prepare_pie_category(y, pie = pie)

    # Get selected categories
    keep <- rownames(ID_Cluster_mat)

    if (length(keep) == 0) {
        stop("no enriched term found...")
    }

    # Prepare similarity matrix
    termsim2 <- fill_termsim(x, keep)

    # Hierarchical clustering
    hc <- hclust(as.dist(1 - termsim2), method = cluster_method)
    nCluster <- min(nCluster, max(1, length(keep)))
    clus <- cutree(hc, nCluster)

    # Prepare data for plotting
    merged_ggData <- merge_compareClusterResult(y)
    rownames(merged_ggData) <- merged_ggData$Description
    d <- data.frame(
        label = names(clus),
        count = merged_ggData[names(clus), "Count"]
    )

    # Create base tree plot
    p <- create_tree_plot(
        hc = hc,
        clus = clus,
        data = d,
        label_format = label_format,
        fontsize_tiplab = fontsize_tiplab,
        fontsize_cladelab = fontsize_cladelab,
        group_color = group_color,
        extend = extend,
        hilight = hilight,
        align = align,
        color_var = color,
        tiplab_offset = tiplab_offset,
        cladelab_offset = cladelab_offset,
        add_tippoint = FALSE # Don't add tip points for compareCluster
    )

    # Add cluster panel based on type
    p <- add_cluster_panel(
        p = p,
        cluster_panel = cluster_panel,
        ID_Cluster_mat = ID_Cluster_mat,
        x = x,
        color = color,
        legend_n = legend_n,
        colnames_angle = colnames_angle,
        hexpand = hexpand
    )

    return(p)
}

#' Add cluster panel to tree plot
#'
#' @param p tree plot
#' @param cluster_panel panel type
#' @param ID_Cluster_mat cluster matrix
#' @param x compareClusterResult object
#' @param color color variable
#' @param legend_n legend items count
#' @param colnames_angle column names angle
#' @param hexpand expand ratio
#' @importFrom rlang sym
#' @noRd
add_cluster_panel <- function(
    p,
    cluster_panel,
    ID_Cluster_mat,
    x,
    color,
    legend_n,
    colnames_angle,
    hexpand
) {
    p_data <- as.data.frame(p$data)
    p_data <- p_data[which(!is.na(p_data$label)), ]
    rownames(p_data) <- p_data$label
    p_data <- p_data[rownames(ID_Cluster_mat), ]

    if (cluster_panel == "pie") {
        # Add pie chart panel
        ID_Cluster_mat$radius <- sqrt(p_data$count / sum(p_data$count))
        ID_Cluster_mat$x <- p_data$x
        ID_Cluster_mat$y <- p_data$y
        ID_Cluster_mat$node <- p_data$node

        p <- p +
            ggnewscale::new_scale_fill() +
            scatterpie::geom_scatterpie(
                aes(x = .data$x, y = .data$y, r = .data$radius),
                data = ID_Cluster_mat,
                cols = colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat) - 4)],
                color = NA
            ) +
            scatterpie::geom_scatterpie_legend(
                ID_Cluster_mat$radius,
                x = 0.8,
                y = 0.1,
                n = legend_n,
                labeller = function(x) round(sum(p_data$count) * x^2)
            ) +
            labs(fill = "Cluster") +
            coord_equal()
    } else if (cluster_panel == "heatMap") {
        # Add heatmap panel
        heatMapData <- as.data.frame(x)
        heatMapData$Cluster <- as.character(heatMapData$Cluster)
        heatMapData <- heatMapData[
            heatMapData$Cluster %in% colnames(ID_Cluster_mat),
        ]
        heatMapData <- heatMapData[
            heatMapData$Description %in% rownames(ID_Cluster_mat),
        ]

        for (i in seq_len(nrow(heatMapData))) {
            ID_Cluster_mat[
                heatMapData[i, "Description"],
                heatMapData[i, "Cluster"]
            ] <- heatMapData[i, color]
        }

        p <- ggtree::gheatmap(
            p + ggnewscale::new_scale_fill(),
            ID_Cluster_mat,
            colnames_angle = colnames_angle,
            width = 0.5
        ) +
            set_enrichplot_color(
                type = "fill",
                transform = "log10",
                name = color
            )
    } else if (cluster_panel == "dotplot") {
        # Add dotplot panel
        require_suggested('ggtreeExtra', 'for `treeplot(cluster_panel = \"dotplot\")`.')
        dotdata <- as.data.frame(x)
        pData <- as.data.frame(p$data)
        paths <- pData$label[order(pData$y, decreasing = TRUE)]
        paths <- paths[!is.na(paths)]
        dotdata <- dotdata[dotdata$Description %in% paths, ]
        dotdata <- dplyr::select(dotdata, "Description", dplyr::everything())

        p <- p +
            ggnewscale::new_scale_colour() +
            ggtreeExtra::geom_fruit(
                data = dotdata,
                geom = geom_point,
                mapping = aes(
                    x = Cluster,
                    y = Description,
                    size = Count,
                    color = .data[[color]]
                ),
                pwidth = 0.06 * ncol(ID_Cluster_mat),
                axis.params = list(
                    axis = "x",
                    text.size = 3,
                    line.alpha = 0,
                    text.angle = colnames_angle
                )
            ) +
            set_enrichplot_color(transform = "log10", name = color)
    }

    return(p + ggtree::hexpand(ratio = hexpand))
}

#' Create tree plot from clustering results
#'
#' @param hc hierarchical clustering result
#' @param clus cluster assignments
#' @param data node data
#' @param label_format label formatting
#' @param fontsize_tiplab tip label font size
#' @param fontsize_cladelab clade label font size
#' @param group_color group colors
#' @param extend extend length
#' @param hilight whether to highlight
#' @param align highlight alignment
#' @param color_var color variable name
#' @param tiplab_offset tip label offset
#' @param cladelab_offset clade label offset
#' @param add_tippoint whether to add tip points (default: TRUE)
#' @noRd
#' @importFrom ggfun %<+%
create_tree_plot <- function(
    hc,
    clus,
    data,
    label_format,
    fontsize_tiplab,
    fontsize_cladelab,
    group_color,
    extend,
    hilight,
    align = 'left',
    color_var,
    size_var = 'size',
    tiplab_offset = 0.2,
    cladelab_offset,
    add_tippoint = TRUE
) {
    # Set colors
    if (is.null(group_color)) {
        require_suggested('scales', 'for `treeplot()`.')
        n_clusters <- length(unique(clus))
        group_color <- scales::hue_pal()(n_clusters)
    }

    # Create base tree
    p <- ggtree(hc, hang = -1, branch.length = "none")

    # Group nodes
    dat <- data.frame(
        name = names(clus),
        cls = paste0("cluster_", as.numeric(clus))
    )
    grp <- apply(table(dat), 2, function(x) names(x[x == 1]))
    clades <- vapply(grp, \(nodes) ggtree::MRCA(p, nodes), numeric(1))
    p <- groupClade(p, clades, "group") +
        aes(color = .data$group) +
        scale_color_manual(
            values = c(group_color, "white"),
            breaks = names(clades)
        )

    # Add tip points and labels
    p <- p %<+% data

    # Add clade labels and highlights
    if (hilight) {
        p <- add_clade_labels(
            p,
            clades,
            label_format,
            fontsize_cladelab,
            group_color,
            extend,
            align,
            offset = cladelab_offset
        )
    }

    if (add_tippoint) {
        p <- p +
            ggnewscale::new_scale_colour() +
            geom_tippoint(aes(
                color = .data[[color_var]],
                size = .data[[size_var]]
            ))
        if (color_var %in% c("pvalue", "qvalue", "p.adjust")) {
            p <- p + set_enrichplot_color(transform = 'log10')
        } else {
            p <- p +
                set_enrichplot_color(
                    colors = rev(get_enrichplot_color(3)),
                )
        }
    }

    p <- p +
        geom_tiplab(
            offset = tiplab_offset,
            hjust = 0,
            size = fontsize_tiplab
        )

    return(p)
}


#' Add clade labels and highlights to tree plot
#'
#' @param p tree plot
#' @param clades clade definitions
#' @param label_format label formatting
#' @param fontsize font size
#' @param group_color group colors
#' @param extend extend length
#' @param align highlight alignment
#' @importFrom ggplot2 scale_fill_manual
#' @noRd
add_clade_labels <- function(
    p,
    clades,
    label_format,
    fontsize,
    group_color,
    extend,
    align,
    offset
) {
    # Prepare clade label data
    df <- data.frame(
        node = as.numeric(clades),
        labels = names(clades),
        cluster = factor(names(clades))
    )

    # Get the tree data to access tip labels
    pdata <- as.data.frame(p$data)
    pdata <- pdata[!is.na(pdata$label), ]

    # Create the required data structure for get_wordcloud
    wordcloud_data <- data.frame(
        name = pdata$label,
        color2 = pdata$group
    )

    # Generate meaningful cluster labels from tip labels
    cluster_labels <- sapply(names(clades), function(cluster_name) {
        get_wordcloud(cluster_name, wordcloud_data, nWords = 4)
    })

    df$labels <- cluster_labels

    # Apply label formatting
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    df$labels <- label_func(df$labels)

    df$color <- group_color

    # Add clade labels and highlights
    p <- p +
        ggnewscale::new_scale_colour() +
        geom_cladelab(
            data = df,
            mapping = aes(
                node = !!sym('node'),
                label = !!sym('labels'),
                color = !!sym('cluster')
            ),
            textcolor = "black",
            extend = extend,
            show.legend = FALSE,
            fontsize = fontsize,
            offset = offset
        ) +
        scale_color_manual(values = group_color, guide = 'none') +
        geom_hilight(
            data = df,
            mapping = aes(node = !!sym('node'), fill = !!sym('cluster')),
            show.legend = FALSE,
            align = align
        ) +
        scale_fill_manual(values = group_color, guide = 'none')

    return(p)
}


#' Fill the upper triangular matrix completely
#'
#' @param x enrichment result
#' @param keep selected categories
#' @return filled similarity matrix
#' @noRd
fill_termsim <- function(x, keep) {
    termsim <- x@termsim[keep, keep]
    termsim[which(is.na(termsim))] <- 0
    termsim2 <- termsim + t(termsim)
    diag(termsim2) <- 1
    return(termsim2)
}
