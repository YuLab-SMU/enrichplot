#' @rdname consensusmap
#' @exportMethod consensusmap
setMethod(
    "consensusmap",
    signature(x = "mnseaResult"),
    function(x, ...) consensusmap_internal(x, ...)
)

#' @rdname consensusmap
#' @exportMethod consensusmap
setMethod(
    "consensusmap",
    signature(x = "nseaResult"),
    function(x, ...) {
        stop("A single nseaResult does not provide multiple networks/conditions. Pass a named list of results to consensusmap().")
    }
)

#' @rdname consensusmap
#' @exportMethod consensusmap
setMethod(
    "consensusmap",
    signature(x = "list"),
    function(x, ...) consensusmap_internal(x, ...)
)

consensusmap_internal <- function(
    x,
    include_rewiring = TRUE,
    ...
) {
    if (is.list(x) && !inherits(x, c("nseaResult", "mnseaResult"))) {
        if (is.null(names(x))) {
            stop("A list input to consensusmap() must be named.")
        }
        contexts <- names(x)
        summaries <- lapply(seq_along(x), function(i) {
            df <- summarize_nsea_mechanism(x[[i]], ...)
            df$context <- contexts[i]
            df
        })
        df <- do.call(rbind, summaries)
    } else if (inherits(x, "mnseaResult")) {
        # Use layer as the context for a single mnseaResult.
        layers <- unique(as.character(x@layer_scores %>% names()))
        summaries <- lapply(layers, function(layer) {
            df <- summarize_nsea_mechanism(x, selected_layer = layer, reference_layer = NULL, ...)
            df$context <- layer
            df
        })
        df <- do.call(rbind, summaries)
    } else {
        stop("x must be a named list or a mnseaResult.")
    }

    if (nrow(df) == 0) {
        stop("No pathways available for consensusmap.")
    }

    df$mechanism_class[is.na(df$mechanism_class)] <- "unknown"
    df$context <- factor(df$context, levels = unique(df$context))
    df$Description <- factor(df$Description, levels = unique(df$Description))

    p <- ggplot(df, aes(x = .data$context, y = .data$Description)) +
        geom_tile(aes(fill = .data$NES)) +
        geom_text(aes(label = .data$mechanism_class), size = 3) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "NES") +
        labs(x = "Context", y = NULL) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    p
}
