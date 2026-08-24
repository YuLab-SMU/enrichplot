#' @rdname mechanismflow
#' @exportMethod mechanismflow
setMethod(
    "mechanismflow",
    signature(x = "mnseaResult"),
    function(x, ...) mechanismflow_internal(x, ...)
)

#' @rdname mechanismflow
#' @exportMethod mechanismflow
setMethod(
    "mechanismflow",
    signature(x = "nseaResult"),
    function(x, ...) mechanismflow_internal(x, ...)
)

#' @rdname mechanismflow
#' @exportMethod mechanismflow
setMethod(
    "mechanismflow",
    signature(x = "list"),
    function(x, ...) mechanismflow_internal(x, ...)
)

mechanismflow_internal <- function(
    x,
    ...
) {
    if (is.list(x) && !inherits(x, c("nseaResult", "mnseaResult"))) {
        if (is.null(names(x))) {
            stop("A list input to mechanismflow() must be named.")
        }
        contexts <- names(x)
        summaries <- lapply(seq_along(x), function(i) {
            df <- summarize_nsea_mechanism(x[[i]], ...)
            df$context <- contexts[i]
            df
        })
        df <- do.call(rbind, summaries)
    } else if (inherits(x, "mnseaResult")) {
        layers <- unique(as.character(x@layer_scores %>% names()))
        summaries <- lapply(layers, function(layer) {
            df <- summarize_nsea_mechanism(x, selected_layer = layer, reference_layer = NULL, ...)
            df$context <- layer
            df
        })
        df <- do.call(rbind, summaries)
    } else if (inherits(x, "nseaResult")) {
        stop("mechanismflow() for a single nseaResult requires a named list of results.")
    } else {
        stop("x must be a named list, nseaResult, or mnseaResult.")
    }

    if (nrow(df) == 0 || length(unique(df$context)) < 2) {
        stop("At least two contexts are required for mechanismflow().")
    }

    df$context <- factor(df$context, levels = unique(df$context))
    df$Description <- factor(df$Description, levels = unique(df$Description))
    df$mechanism_class[is.na(df$mechanism_class)] <- "unknown"
    df$state_num <- as.integer(factor(df$mechanism_class))

    ggplot(df, aes(x = .data$context, y = .data$state_num, group = .data$ID, color = .data$Description)) +
        geom_line(linewidth = 1, alpha = 0.8) +
        geom_point(size = 2) +
        scale_y_continuous(breaks = seq_len(length(unique(df$mechanism_class))), labels = unique(df$mechanism_class)) +
        labs(x = "Context", y = "Mechanism state", color = "Pathway") +
        theme_bw()
}
