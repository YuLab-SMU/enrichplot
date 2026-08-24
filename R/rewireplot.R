#' @rdname rewireplot
#' @exportMethod rewireplot
setMethod(
    "rewireplot",
    signature(x = "nseaResult"),
    function(x, ...) {
        stop("rewireplot() for a single nseaResult requires an explicit reference network/object.")
    }
)

#' @rdname rewireplot
#' @exportMethod rewireplot
setMethod(
    "rewireplot",
    signature(x = "mnseaResult"),
    function(x, ...) rewireplot_internal(x, ...)
)

rewireplot_internal <- function(
    x,
    pathway_id = NULL,
    selected_layer = NULL,
    reference_layer = NULL,
    include_couplings = TRUE,
    node_label = TRUE,
    ...
) {
    if (inherits(x, "mnseaResult") && is.null(reference_layer)) {
        stop("reference_layer must be supplied.")
    }
    if (is.null(pathway_id)) {
        pathway_id <- default_mnsea_pathway_id(x)
    }
    if (is.null(selected_layer)) {
        stop("selected_layer must be supplied.")
    }

    subnet <- fortify_mnsea_subnetwork(
        x,
        pathway_id = pathway_id,
        layer = selected_layer,
        include_couplings = include_couplings
    )
    nodes <- subnet$nodes
    edges <- subnet$edges
    if (nrow(nodes) == 0) {
        stop("No nodes available for the selected pathway/layer.")
    }
    if (nrow(edges) == 0) {
        # Fall back to a simple one-point plot when there are no edges.
        nodes$x <- 0
        nodes$y <- 0
    } else {
        edgelist <- as.data.frame(
            edges[, intersect(c("from", "to", "edge_type", "abs_weight"), colnames(edges)), drop = FALSE],
            stringsAsFactors = FALSE
        )
        node_keys <- intersect(nodes$node_key, unique(c(edgelist$from, edgelist$to)))
        keep_nodes <- nodes$node_key %in% node_keys
        if (sum(keep_nodes) < 1) {
            keep_nodes <- rep(TRUE, nrow(nodes))
        }
        nodes <- nodes[keep_nodes, , drop = FALSE]
        if (nrow(edges) > 0) {
            edges <- edges[edges$from %in% nodes$node_key & edges$to %in% nodes$node_key, , drop = FALSE]
        }
        if (nrow(edges) == 0) {
            nodes$x <- 0
            nodes$y <- 0
        } else {
            g <- igraph::graph_from_data_frame(
                d = edges[, c("from", "to"), drop = FALSE],
                vertices = nodes[, c("node_key", "label", "node_type", "layer", "abs_score", "Feature")],
                directed = FALSE
            )
            lay <- igraph::layout_with_kk(g)
            lay <- as.data.frame(lay)
            colnames(lay) <- c("x", "y")
            lay$node_key <- igraph::V(g)$name
            nodes <- merge(nodes, lay, by = "node_key", all.x = TRUE, sort = FALSE)
        }
    }

    # Attach rewiring status from the helper if reference exists.
    feature_status <- extract_rewiring_features(
        x,
        pathway_id = pathway_id,
        selected_layer = selected_layer,
        reference_layer = reference_layer
    )
    if (nrow(feature_status) > 0) {
        nodes <- merge(
            nodes,
            feature_status,
            by = "Feature",
            all.x = TRUE,
            sort = FALSE
        )
    } else {
        nodes$status <- NA_character_
    }

    p <- ggplot()
    if (nrow(edges) > 0) {
        edge_map <- data.frame(
            x = nodes$x[match(edges$from, nodes$node_key)],
            y = nodes$y[match(edges$from, nodes$node_key)],
            xend = nodes$x[match(edges$to, nodes$node_key)],
            yend = nodes$y[match(edges$to, nodes$node_key)],
            edge_type = edges$edge_type,
            stringsAsFactors = FALSE
        )
        p <- p +
            geom_segment(
                data = edge_map,
                aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, linetype = .data$edge_type),
                color = "grey60",
                linewidth = 0.5
            ) +
            ggplot2::scale_linetype_manual(values = c(membership = "solid", intra = "solid", coupling = "dashed"), guide = "legend")
    }

    node_df <- nodes[nodes$node_type != "pathway", , drop = FALSE]
    p <- p +
        geom_point(
            data = node_df,
            aes(
                x = .data$x,
                y = .data$y,
                size = .data$abs_score,
                color = .data$status,
                shape = .data$node_type
            )
        ) +
        ggplot2::scale_shape_manual(values = c(feature = 19, pathway = 17), guide = "none")
    if (node_label) {
        p <- p + ggrepel::geom_text_repel(
            data = nodes,
            aes(x = .data$x, y = .data$y, label = .data$label),
            inherit.aes = FALSE,
            max.overlaps = 20
        )
    }

    p +
        labs(
            x = NULL,
            y = NULL,
            size = "Feature magnitude",
            color = "Status"
        ) +
        theme_void() +
        theme(legend.position = "right")
}
