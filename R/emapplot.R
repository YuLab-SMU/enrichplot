#' @rdname emapplot
#' @exportMethod emapplot
setMethod(
    "emapplot",
    signature(x = "enrichResult"),
    function(
        x,
        layout = igraph::layout_with_kk,
        coords = NULL,
        showCategory = 30,
        color = "p.adjust",
        size_category = 1,
        min_edge = .2,
        color_edge = "grey",
        size_edge = .5,
        group = NULL,
        group_legend = TRUE,
        node_label = "category",
        node_label_size = 5,
        pie = "equal",
        layer = NULL,
        label_format = 30,
        clusterFunction = stats::kmeans,
        nWords = 4,
        nCluster = NULL,
        ...
    ) {
        emapplot_internal(
            x,
            layout = layout,
            coords = coords,
            showCategory = showCategory,
            color = color,
            size_category = size_category,
            min_edge = min_edge,
            color_edge = color_edge,
            size_edge = size_edge,
            group = group,
            group_legend = group_legend,
            node_label = node_label,
            node_label_size = node_label_size,
            pie = pie,
            layer = layer,
            label_format = label_format,
            clusterFunction = clusterFunction,
            nWords = nWords,
            nCluster = nCluster,
            ...
        )
    }
)

#' @rdname emapplot
#' @exportMethod emapplot
setMethod(
    "emapplot",
    signature(x = "gseaResult"),
    function(
        x,
        layout = igraph::layout_with_kk,
        coords = NULL,
        showCategory = 30,
        color = "p.adjust",
        size_category = 1,
        min_edge = .2,
        color_edge = "grey",
        size_edge = .5,
        group = NULL,
        group_legend = TRUE,
        node_label = "category",
        node_label_size = 5,
        pie = "equal",
        layer = NULL,
        label_format = 30,
        clusterFunction = stats::kmeans,
        nWords = 4,
        nCluster = NULL,
        ...
    ) {
        emapplot_internal(
            x,
            layout = layout,
            coords = coords,
            showCategory = showCategory,
            color = color,
            size_category = size_category,
            min_edge = min_edge,
            color_edge = color_edge,
            size_edge = size_edge,
            group = group,
            group_legend = group_legend,
            node_label = node_label,
            node_label_size = node_label_size,
            pie = pie,
            layer = layer,
            label_format = label_format,
            clusterFunction = clusterFunction,
            nWords = nWords,
            nCluster = nCluster,
            ...
        )
    }
)

#' @rdname emapplot
#' @exportMethod emapplot
setMethod(
    "emapplot",
    signature(x = "compareClusterResult"),
    function(
        x,
        layout = igraph::layout_with_kk,
        coords = NULL,
        showCategory = 30,
        color = "p.adjust",
        size_category = 1,
        min_edge = .2,
        color_edge = "grey",
        size_edge = .5,
        group = NULL,
        group_legend = TRUE,
        node_label = "category",
        node_label_size = 5,
        pie = "equal",
        layer = NULL,
        label_format = 30,
        clusterFunction = stats::kmeans,
        nWords = 4,
        nCluster = NULL,
        ...
    ) {
        emapplot_internal(
            x,
            layout = layout,
            coords = coords,
            showCategory = showCategory,
            color = color,
            size_category = size_category,
            min_edge = min_edge,
            color_edge = color_edge,
            size_edge = size_edge,
            group = group,
            group_legend = group_legend,
            node_label = node_label,
            node_label_size = node_label_size,
            pie = pie,
            layer = layer,
            label_format = label_format,
            clusterFunction = clusterFunction,
            nWords = nWords,
            nCluster = nCluster,
            ...
        )
    }
)

#' @rdname emapplot
#' @exportMethod emapplot
setMethod(
    "emapplot",
    signature(x = "mnseaResult"),
    function(
        x,
        layout = igraph::layout_with_kk,
        coords = NULL,
        showCategory = 30,
        color = "p.adjust",
        size_category = 1,
        min_edge = .2,
        color_edge = "grey",
        size_edge = .5,
        group = NULL,
        group_legend = TRUE,
        node_label = "category",
        node_label_size = 5,
        pie = "equal",
        layer = NULL,
        label_format = 30,
        clusterFunction = stats::kmeans,
        nWords = 4,
        nCluster = NULL,
        ...
    ) {
        emapplot_internal(
            x,
            layout = layout,
            coords = coords,
            showCategory = showCategory,
            color = color,
            size_category = size_category,
            min_edge = min_edge,
            color_edge = color_edge,
            size_edge = size_edge,
            group = group,
            group_legend = group_legend,
            node_label = node_label,
            node_label_size = node_label_size,
            pie = pie,
            layer = layer,
            label_format = label_format,
            clusterFunction = clusterFunction,
            nWords = nWords,
            nCluster = nCluster,
            ...
        )
    }
)


#' @rdname emapplot
#' @param layout igraph layout
#' @param coords optional user-supplied coordinate data.frame with `x`, `y`
#'   and row names matching node labels.
#' @param color Variable used to color enriched terms, e.g. 'pvalue',
#' 'p.adjust' or 'qvalue'.
#' @param size_category relative size of the categories
#' @param min_edge The minimum similarity threshold for whether
#' two nodes are connected, should be between 0 and 1, default value is 0.2.
#' @param color_edge color of the network edge
#' @param size_edge relative size of edge width
#' @param node_label Select which labels to display,
#' one of 'category', 'group', 'all' and 'none'.
#' @param node_label_size size of node label, default is 5.
#' @param pie one of 'equal' or 'Count' to set the slice ratio of the pies
#' @param layer optional layer or layers to retain for `mnseaResult` plots.
#' @param group logical, if TRUE, group the categories.
#' @param group_legend logical, if TRUE, draw a legend for the groups.
#' @param label_format a numeric value sets wrap length, alternatively a custom function to format axis labels.
#' @param clusterFunction clustering method function, such as `stats::kmeans` (default),
#' `cluster::clara`, `cluster::fanny`, or `cluster::pam`.
#' @param nWords Numeric, the number of words in the cluster tags, the default value is 4.
#' @param nCluster Numeric, the number of clusters,
#' the default value is square root of the number of nodes.
#' @importFrom ggplot2 scale_size
#' @importFrom ggtangle geom_edge
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggrepel geom_label_repel
#' @author Guangchuang Yu
prepare_emapplot_data <- function(x, showCategory, color, min_edge, size_edge) {
    ## this path feeds x@termsim straight into the graph builder, so an
    ## unpopulated matrix used to crash with "no 'dimnames' attribute for array"
    has_pairsim(x)
    selected <- select_terms(x, showCategory)
    g <- build_emap_graph(
        enrichDf = selected$result,
        geneSets = selected$geneSets,
        color = color,
        cex_line = size_edge,
        min_edge = min_edge,
        pair_sim = x@termsim
    )
    plot_result <- selected$result
    plot_result$Description <- unname(selected$labels)

    list(
        graph = g,
        geneSet = selected$geneSets,
        result = plot_result
    )
}

prepare_emapplot_mnsea_feature_data <- function(x, ids, layer = NULL) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) == 0) {
        return(data.frame())
    }

    feature_list <- lapply(ids, function(id) {
        fortify_mnsea_contribution(
            x,
            level = "feature",
            pathway_id = id,
            layer = layer
        )
    })
    feature_df <- do.call(rbind, feature_list)
    if (is.null(feature_df) || nrow(feature_df) == 0) {
        return(data.frame())
    }
    rownames(feature_df) <- NULL
    feature_df
}

prepare_mnsea_similarity_data <- function(x, showCategory, layer = NULL) {
    selected <- select_terms(x, showCategory)
    if (nrow(selected$result) == 0) {
        yulab.utils::yulab_abort("No mnsea pathways available for plotting.")
    }

    pathway_df <- fortify_mnsea_contribution(x, level = "pathway", layer = layer)
    feature_df <- prepare_emapplot_mnsea_feature_data(x, ids = selected$ids, layer = layer)
    if (nrow(pathway_df) == 0 || nrow(feature_df) == 0) {
        yulab.utils::yulab_abort("No mnsea layers available after filtering.")
    }

    pathway_df <- pathway_df[pathway_df$ID %in% selected$ids, , drop = FALSE]
    feature_df <- feature_df[feature_df$ID %in% selected$ids, , drop = FALSE]
    retained_ids <- intersect(unique(pathway_df$ID), unique(feature_df$ID))
    if (length(retained_ids) == 0) {
        yulab.utils::yulab_abort("No mnsea pathways available after filtering.")
    }

    keep <- selected$ids %in% retained_ids
    selected$result <- selected$result[keep, , drop = FALSE]
    selected$ids <- selected$ids[keep]
    selected$labels <- selected$labels[keep]
    selected$geneSets <- lapply(selected$ids, function(id) {
        unique(feature_df$Feature[feature_df$ID == id])
    })
    names(selected$geneSets) <- selected$ids
    selected$geneSets <- set_geneSet_labels(selected$geneSets, selected$labels)

    pair_sim <- x@termsim
    method <- x@method
    term_labels <- unname(selected$labels)

    use_cached_termsim <- length(pair_sim) > 0 &&
        !is.null(dim(pair_sim)) &&
        all(term_labels %in% rownames(pair_sim))

    if (use_cached_termsim) {
        pair_sim <- pair_sim[term_labels, term_labels, drop = FALSE]
        if (!nzchar(method)) {
            method <- "JC"
        }
    } else {
        pair_sim <- get_similarity_matrix(
            y = selected$result,
            geneSets = selected$geneSets,
            method = "JC"
        )
        method <- "JC"
    }

    pathway_summary <- stats::aggregate(
        cbind(contribution, share, n_feature) ~ ID + Description,
        data = pathway_df[pathway_df$ID %in% selected$ids, , drop = FALSE],
        FUN = max
    )
    pathway_summary$Description <- as.character(pathway_summary$Description)

    plot_result <- selected$result
    plot_result$Description <- term_labels
    plot_result <- merge(
        plot_result,
        pathway_summary,
        by = c("ID", "Description"),
        all.x = TRUE,
        sort = FALSE
    )

    list(
        result = plot_result,
        geneSet = selected$geneSets,
        pair_sim = pair_sim,
        method = method
    )
}

prepare_emapplot_mnsea_data <- function(x, showCategory, color, min_edge, size_edge, layer = NULL) {
    plot_data <- prepare_mnsea_similarity_data(
        x,
        showCategory = showCategory,
        layer = layer
    )

    g <- build_emap_graph(
        enrichDf = plot_data$result,
        geneSets = plot_data$geneSet,
        color = color,
        cex_line = size_edge,
        min_edge = min_edge,
        pair_sim = plot_data$pair_sim
    )

    list(
        graph = g,
        geneSet = plot_data$geneSet,
        result = plot_data$result
    )
}

emapplot_internal <- function(
    x,
    layout = igraph::layout_with_kk,
    coords = NULL,
    showCategory = 30,
    color = "p.adjust",
    size_category = 1,
    min_edge = .2,
    color_edge = "grey",
    size_edge = .5,
    group = NULL,
    group_legend = TRUE,
    node_label = "category",
    node_label_size = 5,
    pie = "equal",
    layer = NULL,
    label_format = 30,
    clusterFunction = stats::kmeans,
    nWords = 4,
    nCluster = NULL
) {
    if (inherits(x, 'compareClusterResult')) {
        gg <- graph_from_compareClusterResult(
            x,
            showCategory = showCategory,
            color = color,
            min_edge = min_edge,
            size_edge = size_edge
        )
    } else if (inherits(x, 'mnseaResult')) {
        gg <- prepare_emapplot_mnsea_data(
            x,
            showCategory = showCategory,
            color = color,
            min_edge = min_edge,
            size_edge = size_edge,
            layer = layer
        )
    } else {
        gg <- prepare_emapplot_data(x, showCategory, color, min_edge, size_edge)
    }

    g <- gg$graph
    size <- vapply(gg$geneSet, length, FUN.VALUE = numeric(1))
    names(size) <- unname(get_geneSet_labels(gg$geneSet))
    V(g)$size = size[V(g)$name]

    if (!is.null(coords)) {
        coords <- as.data.frame(coords)
        if (!all(c("x", "y") %in% colnames(coords))) {
            yulab.utils::yulab_abort("`coords` must contain `x` and `y` columns.")
        }
        coords <- coords[, c("x", "y"), drop = FALSE]
        if (is.null(rownames(coords))) {
            yulab.utils::yulab_abort("`coords` must use node labels as row names.")
        }
        layout_coords <- coords
        layout <- function(graph) {
            as.matrix(layout_coords[igraph::V(graph)$name, c("x", "y"), drop = FALSE])
        }
    }

    p <- ggplot(g, layout = layout)
    if (igraph::ecount(g) > 0) {
        p <- p + geom_edge(color = color_edge, linewidth = size_edge)
    }

    if (inherits(x, 'compareClusterResult')) {
        p <- add_node_pie(p, gg$data, pie, category_scale = size_category)
    } else {
        if (color %in% names(gg$result)) {
            color_scale <- switch(
                color,
                NES = list(
                    colors = get_enrichplot_color(3),
                    transform = "identity",
                    reverse = TRUE
                ),
                contribution = list(
                    colors = get_enrichplot_color(2),
                    transform = "identity",
                    reverse = TRUE
                ),
                share = list(
                    colors = get_enrichplot_color(2),
                    transform = "identity",
                    reverse = TRUE
                ),
                list(
                    colors = get_enrichplot_color(2),
                    transform = "log10",
                    reverse = TRUE
                )
            )
            p <- p %<+%
                gg$result[, c("Description", color)] +
                geom_point(aes(color = .data[[color]], size = .data$size)) +
                scale_size(
                    range = c(3, 8) * size_category,
                    name = if (inherits(x, "mnseaResult")) "Feature count" else ggplot2::waiver()
                )
            p <- p + set_enrichplot_color(
                colors = color_scale$colors,
                name = if (inherits(x, "mnseaResult")) mnsea_plot_label(color) else color,
                transform = color_scale$transform,
                reverse = color_scale$reverse
            )
            p <- p +
                guides(
                    size = guide_legend(order = 1),
                    color = guide_colorbar(order = 2, reverse = color_scale$reverse)
                )
        } else {
            p <- p %<+%
                gg$result[, "Description", drop = FALSE] +
                geom_point(aes(size = .data$size), color = color) +
                scale_size(
                    range = c(3, 8) * size_category,
                    name = if (inherits(x, "mnseaResult")) "Feature count" else ggplot2::waiver()
                )
        }
    }

    group_enabled <- node_label %in% c("group", "all")
    if (!is.null(group)) {
        group_enabled <- isTRUE(group)
        if (!group_enabled && node_label == "group") {
            node_label <- "none"
        }
    }

    group_label <- FALSE
    if (node_label == "all") {
        group_enabled <- TRUE
        group_label <- TRUE
        node_label <- "category"
    }

    if (group_enabled) {
        if (inherits(x, 'compareClusterResult')) {
            p <- p + ggnewscale::new_scale_fill()
        } #else {
        # p <- p + ggnewscale::new_scale_color()
        #}
        node_data <- groupNode(
            p@data,
            as.data.frame(x),
            nWords,
            clusterFunction = clusterFunction,
            nCluster = nCluster
        )

        p <- p +
            add_ellipse(
                node_data,
                group_legend = group_legend,
                label = group_label
            )
    }

    ## add node label
    if (node_label == "category") {
        p <- p +
            geom_text_repel(
                aes(label = .data$label),
                bg.color = "white",
                bg.r = .1,
                size = node_label_size
            )
    }
    ## add group label
    if (node_label == "group") {
        label_location <- get_label_location(
            node_data = node_data,
            label_format = label_format
        )
        p <- p +
            geom_text_repel(
                aes(x = .data$x, y = .data$y, label = .data$label),
                data = label_location,
                bg.color = "white",
                bg.r = .1,
                size = node_label_size
            )
    }

    p +
        coord_equal() +
        guides(
            size = guide_legend(order = 1),
            color = guide_colorbar(order = 2)
        )
}

graph_from_compareClusterResult <- function(
    x,
    showCategory = 30,
    color = "p.adjust",
    min_edge = .2,
    size_edge = .5
) {
    ## x@termsim goes straight into the graph builder below, so an unpopulated
    ## matrix used to crash with "no 'dimnames' attribute for array"
    has_pairsim(x)
    d <- tidy_compareCluster(x, showCategory)
    mergedEnrichDf <- merge_compareClusterResult(d)
    term_labels <- unname(get_term_labels(mergedEnrichDf, mergedEnrichDf$ID))
    label_map <- stats::setNames(term_labels, mergedEnrichDf$ID)
    d$Description <- unname(label_map[d$ID])
    mergedEnrichDf$Description <- term_labels
    gs <- setNames(
        strsplit(as.character(mergedEnrichDf$geneID), "/", fixed = TRUE),
        mergedEnrichDf$ID
    )

    g <- build_emap_graph(
        enrichDf = mergedEnrichDf,
        geneSets = gs,
        color = color,
        cex_line = size_edge,
        min_edge = min_edge,
        pair_sim = x@termsim
    )
    return(list(graph = g, geneSet = gs, data = d))
}
