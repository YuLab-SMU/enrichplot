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
    fill_var = c("NES", "delta_NES"),
    size_var = c("rewiring_score", "leading_edge_overlap"),
    include_rewiring = TRUE,
    label = c("mechanism_class", "rewiring_score", "none"),
    reference = NULL,
    ...
) {
    fill_var <- match.arg(fill_var)
    size_var <- match.arg(size_var)
    label <- match.arg(label)

    if (is.list(x) && !inherits(x, c("nseaResult", "mnseaResult"))) {
        if (is.null(names(x))) {
            stop("A list input to consensusmap() must be named.")
        }
        contexts <- names(x)
        if (is.null(reference)) {
            reference <- x[[1]]
        }
        summaries <- lapply(seq_along(x), function(i) {
            df <- summarize_nsea_mechanism(
                x[[i]],
                reference = reference,
                ...
            )
            df$context <- contexts[i]
            df
        })
        df <- do.call(rbind, summaries)
    } else if (inherits(x, "mnseaResult")) {
        # Use layer as the context for a single mnseaResult.
        layers <- unique(as.character(x@layer_scores %>% names()))
        summaries <- lapply(layers, function(layer) {
            df <- summarize_nsea_mechanism(
                x,
                reference = reference,
                selected_layer = layer,
                reference_layer = NULL,
                ...
            )
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

    p <- ggplot(df, aes(x = .data$context, y = .data$Description))

    if (fill_var %in% colnames(df) && !all(is.na(df[[fill_var]]))) {
        p <- p +
            geom_tile(aes(fill = .data[[fill_var]])) +
            scale_fill_gradient2(
                low = "blue",
                mid = "white",
                high = "red",
                name = if (fill_var == "delta_NES") "Delta NES" else "NES"
            )
    } else {
        p <- p + geom_tile(fill = "grey90")
    }

    if (include_rewiring && size_var %in% colnames(df) && !all(is.na(df[[size_var]]))) {
        p <- p +
            geom_point(
                aes(
                    size = .data[[size_var]],
                    color = .data$mechanism_class
                ),
                shape = 21,
                fill = "white",
                alpha = 0.8
            ) +
            scale_size_continuous(
                range = c(2, 8),
                name = if (size_var == "rewiring_score") "Rewiring score" else "Leading-edge overlap"
            )
    }

    if (label != "none") {
        if (label == "mechanism_class") {
            p <- p + geom_text(aes(label = .data$mechanism_class), size = 3)
        } else {
            p <- p + geom_text(
                aes(label = ifelse(
                    is.na(.data[[size_var]]),
                    "NA",
                    sprintf("%.2f", .data[[size_var]])
                )),
                size = 3
            )
        }
    }

    p +
        labs(
            x = "Context",
            y = NULL,
            color = "Mechanism"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
}
