#' @rdname phaseplot
#' @exportMethod phaseplot
setMethod(
    "phaseplot",
    signature(x = "nseaResult"),
    function(
        x,
        reference = NULL,
        selected_layer = NULL,
        reference_layer = NULL,
        x_axis = NULL,
        size_var = c("leading_edge_size", "leading_edge_overlap"),
        showCategory = 30,
        ...
    ) {
        phaseplot_internal(
            x,
            reference = reference,
            selected_layer = selected_layer,
            reference_layer = reference_layer,
            x_axis = x_axis,
            size_var = size_var,
            showCategory = showCategory,
            ...
        )
    }
)

#' @rdname phaseplot
#' @exportMethod phaseplot
setMethod(
    "phaseplot",
    signature(x = "mnseaResult"),
    function(
        x,
        reference = NULL,
        selected_layer = NULL,
        reference_layer = NULL,
        x_axis = NULL,
        size_var = c("leading_edge_size", "leading_edge_overlap"),
        showCategory = 30,
        ...
    ) {
        phaseplot_internal(
            x,
            reference = reference,
            selected_layer = selected_layer,
            reference_layer = reference_layer,
            x_axis = x_axis,
            size_var = size_var,
            showCategory = showCategory,
            ...
        )
    }
)

phaseplot_internal <- function(
    x,
    reference = NULL,
    selected_layer = NULL,
    reference_layer = NULL,
    x_axis = NULL,
    size_var = c("leading_edge_size", "leading_edge_overlap"),
    showCategory = 30,
    ...
) {
    x_axis_specified <- !is.null(x_axis)
    if (!is.null(x_axis)) {
        x_axis <- match.arg(x_axis, c("delta_NES", "NES"))
    } else {
        x_axis <- "NES"
    }
    size_var <- match.arg(size_var)

    summary_df <- summarize_nsea_mechanism(
        x,
        reference = reference,
        selected_layer = selected_layer,
        reference_layer = reference_layer
    )
    if (nrow(summary_df) == 0) {
        stop("No mnsea pathways available for phaseplot.")
    }

    if (is.numeric(showCategory) && nrow(summary_df) > showCategory) {
        summary_df <- summary_df[seq_len(showCategory), , drop = FALSE]
    } else if (!is.numeric(showCategory)) {
        keep <- summary_df$ID %in% showCategory | summary_df$Description %in% showCategory
        summary_df <- summary_df[keep, , drop = FALSE]
    }

    if (x_axis == "NES" && !x_axis_specified && !is.null(reference) && !all(is.na(summary_df$delta_NES))) {
        x_axis <- "delta_NES"
    }
    x_col <- if (x_axis == "delta_NES") "delta_NES" else "NES"
    if (x_axis == "delta_NES" && all(is.na(summary_df$delta_NES))) {
        stop("No `delta_NES` available for phaseplot(). Supply `reference` or use `x_axis = \"NES\"`.")
    }
    if (x_axis == "NES" && all(is.na(summary_df$NES))) {
        stop("No `NES` available for phaseplot().")
    }

    p <- ggplot(
        summary_df,
        aes(
            x = .data[[x_col]],
            y = .data$rewiring_score,
            size = .data[[size_var]],
            color = .data$mechanism_class
        )
    ) +
        geom_point(alpha = 0.8) +
        labs(
            x = if (x_axis == "delta_NES") "Enrichment shift (delta NES)" else "Enrichment strength (NES)",
            y = "Rewiring score (1 - leading-edge overlap)",
            size = if (size_var == "leading_edge_size") "Leading edge size" else "Leading edge overlap",
            color = "Mechanism"
        ) +
        theme_bw()

    if (!any(!is.na(summary_df$mechanism_class))) {
        p <- p + scale_color_manual(values = "grey50", na.value = "grey50")
    }

    p
}
