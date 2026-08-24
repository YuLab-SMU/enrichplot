#' @rdname phaseplot
#' @exportMethod phaseplot
setMethod(
    "phaseplot",
    signature(x = "nseaResult"),
    function(x, ...) phaseplot_internal(x, ...)
)

#' @rdname phaseplot
#' @exportMethod phaseplot
setMethod(
    "phaseplot",
    signature(x = "mnseaResult"),
    function(x, ...) phaseplot_internal(x, ...)
)

phaseplot_internal <- function(
    x,
    selected_layer = NULL,
    reference_layer = NULL,
    showCategory = 30,
    ...
) {
    summary_df <- summarize_nsea_mechanism(
        x,
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

    p <- ggplot(
        summary_df,
        aes(
            x = .data$NES,
            y = .data$rewiring_score,
            size = .data$leading_edge_size,
            color = .data$mechanism_class
        )
    ) +
        geom_point(alpha = 0.8) +
        labs(
            x = "Enrichment strength (NES)",
            y = "Rewiring score (1 - leading-edge overlap)",
            size = "Leading edge size",
            color = "Mechanism"
        ) +
        theme_bw()

    if (!any(!is.na(summary_df$mechanism_class))) {
        p <- p + scale_color_manual(values = "grey50", na.value = "grey50")
    }

    p
}
