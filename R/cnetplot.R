#' Category-Gene-Network Plot
#'
#' Category-gene-network plot
#' @rdname cnetplot
#' @param x input object
#' @param layout network layout
#' @param showCategory number of categories to display or a vector of terms.
#' @param color_category color of category nodes
#' @param size_category relative size of the category nodes
#' @param color_item color of item nodes
#' @param size_item relative size of the item nodes (e.g., genes)
#' @param color_edge color of edge
#' @param size_edge relative size of edge
#' @param categorySizeBy An expression (e.g., `itemNum`, `p.adjust`) or a formula
#'   (e.g., `~ -log10(p.adjust)`) to set the category node size. For
#'   `compareClusterResult`, this controls the category pie size.
#' @param node_label one of 'all', 'none', 'category', 'item', 'exclusive' or 'share'.
#' 'exclusive' labels genes that uniquely belong to categories; 'share' labels genes that are shared between categories.
#' @param foldChange numeric values to color the item (e.g., fold change of gene expression values)
#' @param fc_threshold threshold for filtering genes by absolute fold change (e.g., fc_threshold = 1 keeps only genes with |foldChange| > 1).
#' @param hilight selected categories to be highlighted
#' @param hilight_alpha transparency value for non-highlighted items
#' @param split apply `showCategory` to each category specified by `split` for
#'   `compareClusterResult`, e.g. `ONTOLOGY`, `category` or `intersect`.
#' @param includeAll logical value passed to `fortify()` when selecting terms
#'   from a `compareClusterResult`.
#' @param pathway_id optional pathway ID for `mnseaResult` subnetworks.
#' @param layer optional layer or layers to retain for `mnseaResult` plots.
#' @param include_couplings logical, whether inter-layer coupling edges should
#'   be kept in `mnseaResult` network plots.
#' @param include_isolated logical, whether isolated feature nodes should be
#'   kept in `mnseaResult` subnetworks.
#' @param ... additional parameters
#' @importFrom ggtangle cnetplot
#' @seealso
#' [cnetplot][ggtangle::cnetplot]
prepare_cnetplot_data <- function(x, showCategory, foldChange) {
    selected <- select_terms(x, showCategory)
    geneSets <- attach_result_attributes(selected$geneSets, selected$result)

    list(
        geneSets = geneSets,
        plot_geneSets = as_plot_geneSets(geneSets),
        foldChange = fc_readable(x, foldChange)
    )
}

#' @rdname cnetplot
#' @method cnetplot enrichResult
#' @export
cnetplot.enrichResult <- function(
    x,
    layout = igraph::layout_with_kk,
    showCategory = 5,
    color_category = "#E5C494",
    size_category = 1,
    color_item = "#B3B3B3",
    size_item = 1,
    color_edge = "grey",
    size_edge = .5,
    categorySizeBy = ~itemNum,
    node_label = "all",
    foldChange = NULL,
    fc_threshold = NULL,
    hilight = "none",
    hilight_alpha = .3,
    ...
) {
    plot_data <- prepare_cnetplot_data(x, showCategory, foldChange)

    args <- list(...)
    plot_args <- list(
        x = plot_data$plot_geneSets,
        layout = layout,
        showCategory = showCategory,
        foldChange = plot_data$foldChange,
        fc_threshold = fc_threshold,
        color_category = color_category,
        size_category = size_category,
        color_item = color_item,
        size_item = size_item,
        color_edge = color_edge,
        size_edge = size_edge,
        node_label = node_label,
        hilight = hilight,
        hilight_alpha = hilight_alpha,
        categorySizeBy = categorySizeBy
    )
    
    final_args <- c(plot_args, args)
    
    p <- do.call(cnetplot, final_args)

    p <- p +
        set_enrichplot_color(
            colors = get_enrichplot_color(3),
            name = "fold change",
            transform = 'identity'
        )
    if (!is.null(plot_data$foldChange)) {
        p <- p +
            guides(
                size = guide_legend(order = 1),
                color = guide_colorbar(order = 2)
            )
    }

    return(p + guides(alpha = "none"))
}

#' @rdname cnetplot
#' @method cnetplot gseaResult
#' @export
cnetplot.gseaResult <- cnetplot.enrichResult

#' @rdname cnetplot
#' @details For `mnseaResult`, `cnetplot()` visualizes a pathway-specific
#' multilayer subnetwork returned by `fortify_mnsea_subnetwork()`. Feature nodes
#' are grouped by layer, a synthetic pathway node is added as the anchor, and
#' optional coupling edges connect features across layers.
prepare_mnsea_cnetplot_data <- function(
    x,
    pathway_id,
    layer,
    include_couplings,
    include_isolated = FALSE
) {
    subnet <- fortify_mnsea_subnetwork(
        x,
        pathway_id = pathway_id,
        layer = layer,
        include_couplings = include_couplings,
        include_isolated = include_isolated
    )
    nodes <- subnet$nodes
    edges <- subnet$edges

    if (nrow(nodes) == 0) {
        yulab.utils::yulab_abort("No mnsea subnetwork nodes available for plotting.")
    }

    if (!"abs_score" %in% colnames(nodes)) {
        nodes$abs_score <- 1
    }
    feature_max <- max(nodes$abs_score[nodes$node_type == "feature"], na.rm = TRUE)
    if (!is.finite(feature_max) || feature_max <= 0) {
        feature_max <- 1
    }
    nodes$plot_size <- ifelse(
        nodes$node_type == "pathway",
        feature_max,
        nodes$abs_score
    )
    nodes$plot_size[!is.finite(nodes$plot_size) | nodes$plot_size <= 0] <- 1

    if (nrow(edges) == 0) {
        vertices <- nodes
        vertices$name <- vertices$node_key
        vertices <- vertices[, c("name", setdiff(colnames(vertices), "name")), drop = FALSE]
        g <- igraph::graph_from_data_frame(
            data.frame(from = character(0), to = character(0)),
            directed = FALSE,
            vertices = vertices
        )
        return(list(graph = g, nodes = nodes, edges = edges, pathway = subnet$pathway))
    }

    vertices <- nodes
    vertices$name <- vertices$node_key
    vertices <- vertices[, c("name", setdiff(colnames(vertices), "name")), drop = FALSE]
    g <- igraph::graph_from_data_frame(
        edges[, c("from", "to"), drop = FALSE],
        directed = FALSE,
        vertices = vertices
    )
    edge_match <- match(
        paste(igraph::ends(g, igraph::E(g))[, 1], igraph::ends(g, igraph::E(g))[, 2]),
        paste(edges$from, edges$to)
    )
    if (anyNA(edge_match)) {
        reverse_match <- match(
            paste(igraph::ends(g, igraph::E(g))[, 2], igraph::ends(g, igraph::E(g))[, 1]),
            paste(edges$from, edges$to)
        )
        edge_match[is.na(edge_match)] <- reverse_match[is.na(edge_match)]
    }
    igraph::E(g)$edge_type <- edges$edge_type[edge_match]
    igraph::E(g)$abs_weight <- edges$abs_weight[edge_match]

    list(graph = g, nodes = nodes, edges = edges, pathway = subnet$pathway)
}

normalize_mnsea_edge_width <- function(x, size_edge) {
    if (length(x) == 0) {
        return(numeric(0))
    }
    x <- as.numeric(x)
    xmin <- suppressWarnings(min(x, na.rm = TRUE))
    xmax <- suppressWarnings(max(x, na.rm = TRUE))

    if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin) {
        return(rep(size_edge, length(x)))
    }

    size_edge * (0.8 + (x - xmin) / (xmax - xmin) * 1.2)
}

select_mnsea_feature_labels <- function(
    feature_nodes,
    membership = NULL,
    max_labels = 8L,
    prefer_share = TRUE
) {
    if (nrow(feature_nodes) == 0) {
        return(feature_nodes)
    }

    label_nodes <- feature_nodes
    if (!is.null(membership)) {
        label_nodes <- label_nodes[
            label_nodes$membership_class %in% membership,
            ,
            drop = FALSE
        ]
    }
    if (nrow(label_nodes) == 0) {
        return(label_nodes)
    }

    if (prefer_share) {
        membership_rank <- c(share = 0L, exclusive = 1L)
        label_nodes$membership_rank <- membership_rank[label_nodes$membership_class]
        label_nodes$membership_rank[is.na(label_nodes$membership_rank)] <- 2L
    } else {
        label_nodes$membership_rank <- 0L
    }

    label_nodes <- label_nodes[
        order(label_nodes$membership_rank, -label_nodes$plot_size, label_nodes$Feature),
        ,
        drop = FALSE
    ]
    label_nodes <- label_nodes[!duplicated(label_nodes$Feature), , drop = FALSE]

    if (!is.null(max_labels) && nrow(label_nodes) > max_labels) {
        label_nodes <- label_nodes[seq_len(max_labels), , drop = FALSE]
    }

    label_nodes$membership_rank <- NULL
    label_nodes
}

select_mnsea_label_data <- function(node_label, node_data, pathway_nodes, feature_nodes) {
    feature_label_nodes <- select_mnsea_feature_labels(feature_nodes, max_labels = 8L)
    feature_label_nodes_quiet <- select_mnsea_feature_labels(feature_nodes, max_labels = 6L)
    share_label_nodes <- select_mnsea_feature_labels(
        feature_nodes,
        membership = "share",
        max_labels = 6L
    )
    exclusive_label_nodes <- select_mnsea_feature_labels(
        feature_nodes,
        membership = "exclusive",
        max_labels = 6L
    )

    switch(
        node_label,
        category = list(pathway = pathway_nodes, feature = feature_nodes[0, , drop = FALSE]),
        item = list(pathway = pathway_nodes[0, , drop = FALSE], feature = feature_label_nodes),
        exclusive = list(pathway = pathway_nodes[0, , drop = FALSE], feature = exclusive_label_nodes),
        share = list(pathway = pathway_nodes[0, , drop = FALSE], feature = share_label_nodes),
        all = list(pathway = pathway_nodes, feature = feature_label_nodes_quiet),
        list(pathway = pathway_nodes, feature = feature_label_nodes_quiet)
    )
}

add_mnsea_label_layers <- function(p, label_data) {
    if (nrow(label_data$pathway) > 0) {
        p <- p +
            geom_text_repel(
                data = label_data$pathway,
                aes(x = .data$x, y = .data$y, label = .data$label),
                seed = 1,
                size = 4.2,
                fontface = "bold",
                box.padding = 0.5,
                point.padding = 0.35,
                min.segment.length = 0,
                segment.color = "grey60",
                bg.color = "white",
                bg.r = 0.12
            )
    }

    if (nrow(label_data$feature) > 0) {
        p <- p +
            geom_text_repel(
                data = label_data$feature,
                aes(x = .data$x, y = .data$y, label = .data$label),
                seed = 1,
                size = 3.2,
                box.padding = 0.25,
                point.padding = 0.15,
                min.segment.length = 0,
                segment.color = "grey75",
                bg.color = "white",
                bg.r = 0.08
            )
    }

    p
}

#' @rdname cnetplot
#' @method cnetplot mnseaResult
#' @export
cnetplot.mnseaResult <- function(
    x,
    layout = igraph::layout_with_kk,
    showCategory = 5,
    color_category = "#E5C494",
    size_category = 1,
    color_item = "#B3B3B3",
    size_item = 1,
    color_edge = "grey",
    size_edge = .5,
    categorySizeBy = ~itemNum,
    node_label = "all",
    foldChange = NULL,
    fc_threshold = NULL,
    hilight = "none",
    hilight_alpha = .3,
    pathway_id = NULL,
    layer = NULL,
    include_couplings = TRUE,
    ...
) {
    plot_data <- prepare_mnsea_cnetplot_data(
        x = x,
        pathway_id = pathway_id,
        layer = layer,
        include_couplings = include_couplings
    )
    g <- plot_data$graph
    edge_type_labels <- c(
        membership = "Pathway membership",
        intra = "Within-layer",
        coupling = "Cross-layer coupling"
    )
    sign_colors <- c(
        activated = "#D73027",
        suppressed = "#4575B4",
        neutral = "grey40"
    )
    edge_types <- as.character(igraph::E(g)$edge_type)
    igraph::E(g)$edge_type_label <- factor(
        edge_type_labels[edge_types],
        levels = unname(edge_type_labels)
    )
    igraph::E(g)$edge_width <- normalize_mnsea_edge_width(
        igraph::E(g)$abs_weight,
        size_edge = size_edge
    )

    p <- ggplot(g, layout = layout) +
        geom_edge(
            aes(
                linewidth = .data$edge_width,
                linetype = .data$edge_type_label
            ),
            color = color_edge,
            alpha = 0.7
        )

    node_data <- p$data
    node_data$membership_class <- NA_character_
    node_type_labels <- c(
        feature = "Feature node",
        pathway = "Pathway node"
    )
    node_data$node_type_label <- factor(
        node_type_labels[node_data$node_type],
        levels = unname(node_type_labels)
    )
    feature_nodes <- node_data[node_data$node_type == "feature", , drop = FALSE]
    pathway_nodes <- node_data[node_data$node_type == "pathway", , drop = FALSE]
    if (nrow(feature_nodes) > 0) {
        feature_freq <- table(feature_nodes$Feature)
        feature_nodes$membership_class <- ifelse(
            feature_freq[feature_nodes$Feature] > 1,
            "share",
            "exclusive"
        )
    } else {
        feature_nodes$membership_class <- character(0)
    }
    node_data[node_data$node_type == "feature", "membership_class"] <- feature_nodes$membership_class
    pathway_nodes <- node_data[node_data$node_type == "pathway", , drop = FALSE]
    feature_nodes <- node_data[node_data$node_type == "feature", , drop = FALSE]
    layer_levels <- unique(plot_data$nodes$layer[plot_data$nodes$node_type == "feature"])
    layer_levels <- layer_levels[!is.na(layer_levels)]
    if (length(layer_levels) > 0) {
        feature_nodes$layer <- factor(feature_nodes$layer, levels = layer_levels)
    }

    p <- p +
        geom_point(
            data = feature_nodes,
            aes(
                x = .data$x,
                y = .data$y,
                size = .data$plot_size,
                fill = .data$layer,
                color = .data$sign,
                shape = .data$node_type_label
            ),
            stroke = 1
        ) +
        geom_point(
            data = pathway_nodes,
            aes(x = .data$x, y = .data$y, shape = .data$node_type_label),
            size = 6 * size_category,
            fill = color_category,
            color = "black",
            stroke = 1.2
        ) +
        ggplot2::scale_size_continuous(
            range = c(3, 8),
            name = "Feature magnitude"
        ) +
        ggplot2::scale_linewidth_identity() +
        ggplot2::scale_shape_manual(
            values = c(
                "Feature node" = 21,
                "Pathway node" = 23
            ),
            name = "Node type",
            drop = FALSE
        ) +
        ggplot2::scale_color_manual(
            values = sign_colors,
            name = "Feature sign",
            drop = FALSE
        ) +
        ggplot2::scale_fill_discrete(name = "Layer", drop = FALSE) +
        ggplot2::scale_linetype_manual(
            values = c(
                "Pathway membership" = "solid",
                "Within-layer" = "solid",
                "Cross-layer coupling" = "dashed"
            ),
            drop = FALSE,
            name = "Edge type"
        ) +
        guides(
            linewidth = "none",
            linetype = guide_legend(order = 1),
            shape = guide_legend(order = 2),
            fill = guide_legend(order = 3),
            color = guide_legend(order = 4),
            size = guide_legend(order = 5)
        )

    if (node_label != "none") {
        label_data <- select_mnsea_label_data(
            node_label = node_label,
            node_data = node_data,
            pathway_nodes = pathway_nodes,
            feature_nodes = feature_nodes
        )
        p <- add_mnsea_label_layers(p, label_data)
    }

    p
}

#' @rdname cnetplot
#' @param pie one of 'equal' or 'Count' to set the slice ratio of the pies
#' @method cnetplot compareClusterResult
#' @export
cnetplot.compareClusterResult <- function(
    x,
    layout = igraph::layout_with_kk,
    showCategory = 5,
    color_category = "#E5C494",
    size_category = 1,
    color_item = "#B3B3B3",
    size_item = 1,
    color_edge = "grey",
    size_edge = .5,
    categorySizeBy = ~itemNum,
    node_label = "all",
    foldChange = NULL,
    fc_threshold = NULL,
    hilight = "none",
    hilight_alpha = .3,
    pie = "equal",
    split = NULL,
    includeAll = TRUE,
    ...
) {
    category_size_quo <- rlang::enquo(categorySizeBy)
    d <- tidy_compareCluster(
        x,
        showCategory = showCategory,
        split = split,
        includeAll = includeAll
    )
    d$Description <- unname(get_term_labels(x, d$ID))
    y <- split(d$geneID, d$Description)
    gs <- lapply(y, function(item) unique(unlist(strsplit(item, split = "/"))))
    category_size <- compute_comparecluster_category_size(d, category_size_quo)

    p <- cnetplot(
        gs,
        layout = layout,
        showCategory = names(gs),
        foldChange = foldChange,
        fc_threshold = fc_threshold,
        color_category = color_category,
        size_category = 0,
        color_item = color_item,
        size_item = 0,
        color_edge = color_edge,
        size_edge = size_edge,
        node_label = "none",
        hilight = hilight,
        hilight_alpha = hilight_alpha,
        ...
    )

    p <- add_node_pie(
        p,
        d,
        pie,
        category_scale = size_category,
        item_scale = size_item,
        category_size = category_size
    )

    p <- p + geom_cnet_label(node_label = node_label)

    return(p)
}


#' @importFrom ggplot2 coord_fixed
add_node_pie <- function(
    p,
    d,
    pie = "equal",
    category_scale = 1,
    item_scale = 1,
    category_size = NULL
) {
    require_suggested('tidyr', 'for `cnetplot()`.')

    ## category nodes
    dd <- d[, c('Cluster', 'Description', 'Count')]
    if (nrow(dd) > 0) {
        dd <- stats::aggregate(
            Count ~ Cluster + Description,
            data = dd,
            FUN = sum
        )
    }
    default_size <- sapply(split(dd$Count, dd$Description), sum)
    if (is.null(category_size)) {
        category_size <- default_size
    } else {
        category_size <- category_size[names(default_size)]
    }
    if (pie == "equal") {
        dd$Count <- 1
    }
    dd <- tidyr::pivot_wider(
        dd,
        names_from = "Cluster",
        values_from = "Count",
        values_fill = list(Count = 0),
        values_fn = sum
    )
    normalized_category_size <- normalize_comparecluster_radius(category_size)
    dd$pathway_radius <- normalized_category_size[dd$Description] * category_scale

    ## gene nodes
    y <- split(d$geneID, d$Cluster)
    gs <- lapply(y, function(item) unique(unlist(strsplit(item, split = "/"))))
    dg <- yulab.utils::ls2df(gs) |> setNames(c("Cluster", "Description")) # second column is geneID
    dg$Count <- 1
    dg <- tidyr::pivot_wider(
        dg,
        names_from = "Cluster",
        values_from = "Count",
        values_fill = list(Count = 0),
        values_fn = sum
    )
    dg$pathway_radius <- .05 * item_scale

    d2 <- rbind(dd, dg)
    node_layout <- p$data
    node_label_col <- if ("name" %in% colnames(node_layout)) {
        "name"
    } else {
        "label"
    }
    node_layout <- data.frame(
        Description = as.character(node_layout[[node_label_col]]),
        x = node_layout$x,
        y = node_layout$y,
        stringsAsFactors = FALSE
    )
    d2 <- merge(
        d2,
        node_layout,
        by = "Description",
        all.x = FALSE,
        all.y = FALSE,
        sort = FALSE
    )

    p <- p +
        scatterpie::geom_scatterpie(
            data = d2,
            aes(
                x = .data$x,
                y = .data$y,
                r = .data$pathway_radius
            ),
            cols = as.character(unique(d$Cluster)),
            legend_name = "Cluster",
            color = NA,
            inherit.aes = FALSE
        ) +
        coord_fixed() +
        guides(size = "none")

    if (nrow(d2) > 0 && any(dd$pathway_radius > 0)) {
        p <- p + scatterpie::geom_scatterpie_legend(
            unique(dd$pathway_radius),
            x = min(d2$x),
            y = min(d2$y),
            n = 3,
            labeller = function(x) {
                format(signif(x / category_scale * sum(category_size), 3), trim = TRUE)
            }
        )
    }

    return(p)
}

compute_comparecluster_category_size <- function(d, categorySizeBy) {
    term_names <- unique(as.character(d$Description))
    gene_sets <- split(d$geneID, d$Description)
    item_num <- vapply(
        gene_sets[term_names],
        function(item) {
            length(unique(unlist(strsplit(item, split = "/"))))
        },
        FUN.VALUE = numeric(1)
    )

    term_df <- data.frame(
        Description = term_names,
        itemNum = item_num,
        Count = vapply(
            split(d$Count, d$Description)[term_names],
            sum,
            FUN.VALUE = numeric(1)
        ),
        stringsAsFactors = FALSE
    )

    numeric_cols <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], "Count")
    if (length(numeric_cols) > 0) {
        for (col in numeric_cols) {
            term_df[[col]] <- summarize_comparecluster_numeric_column(
                d[[col]],
                d$Description,
                term_names,
                col
            )
        }
    }

    if (inherits(categorySizeBy, "quosure")) {
        category_size_expr <- rlang::quo_get_expr(categorySizeBy)
        category_size_env <- rlang::quo_get_env(categorySizeBy)
    } else {
        category_size_expr <- substitute(categorySizeBy)
        category_size_env <- parent.frame()
    }
    if (rlang::is_formula(category_size_expr)) {
        category_size_expr <- rlang::f_rhs(category_size_expr)
    }

    category_size <- rlang::eval_tidy(
        rlang::new_quosure(category_size_expr, env = category_size_env),
        data = term_df
    )
    if (!is.numeric(category_size)) {
        stop("`categorySizeBy` must evaluate to a numeric vector.")
    }
    if (length(category_size) == 1) {
        category_size <- rep(category_size, nrow(term_df))
    }
    if (length(category_size) != nrow(term_df)) {
        stop("`categorySizeBy` must return a scalar or one value per category.")
    }
    if (anyNA(category_size) || any(!is.finite(category_size))) {
        stop("`categorySizeBy` returned non-finite values.")
    }
    if (any(category_size < 0)) {
        stop("`categorySizeBy` must be non-negative for pie radius scaling.")
    }
    if (sum(category_size) <= 0) {
        stop("`categorySizeBy` must produce at least one positive value.")
    }

    stats::setNames(category_size, term_df$Description)
}

summarize_comparecluster_numeric_column <- function(values, groups, group_names, column_name) {
    split_values <- split(values, groups)
    vapply(
        split_values[group_names],
        function(x) {
            x <- x[is.finite(x)]
            if (length(x) == 0) {
                return(NA_real_)
            }
            if (column_name %in% c("pvalue", "p.adjust", "qvalue")) {
                return(min(x))
            }
            unique_values <- unique(x)
            if (length(unique_values) == 1) {
                return(unique_values)
            }
            NA_real_
        },
        FUN.VALUE = numeric(1)
    )
}

normalize_comparecluster_radius <- function(x) {
    total <- sum(x)
    if (!is.finite(total) || total <= 0) {
        stop("Category pie size scaling requires a positive total size.")
    }
    x / total
}


tidy_compareCluster <- function(
    x,
    showCategory,
    split = NULL,
    includeAll = TRUE
) {
    d <- fortify(
        x,
        showCategory = showCategory,
        includeAll = includeAll,
        split = split
    )
    d$Cluster <- sub("\n.*", "", d$Cluster)

    if ("core_enrichment" %in% colnames(d)) {
        ## for GSEA result
        d$geneID <- d$core_enrichment
    }
    return(d)
}
