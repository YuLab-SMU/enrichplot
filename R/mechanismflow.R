#' @rdname mechanismflow
#' @exportMethod mechanismflow
setMethod(
    "mechanismflow",
    signature(x = "mnseaResult"),
    function(
        x,
        reference = NULL,
        flow_var = c("NES", "delta_NES", "leading_edge_size"),
        ...
    ) {
        mechanismflow_internal(
            x,
            reference = reference,
            flow_var = flow_var,
            ...
        )
    }
)

#' @rdname mechanismflow
#' @exportMethod mechanismflow
setMethod(
    "mechanismflow",
    signature(x = "nseaResult"),
    function(
        x,
        reference = NULL,
        flow_var = c("NES", "delta_NES", "leading_edge_size"),
        ...
    ) {
        mechanismflow_internal(
            x,
            reference = reference,
            flow_var = flow_var,
            ...
        )
    }
)

#' @rdname mechanismflow
#' @exportMethod mechanismflow
setMethod(
    "mechanismflow",
    signature(x = "list"),
    function(
        x,
        reference = NULL,
        flow_var = c("NES", "delta_NES", "leading_edge_size"),
        ...
    ) {
        mechanismflow_internal(
            x,
            reference = reference,
            flow_var = flow_var,
            ...
        )
    }
)

mechanismflow_internal <- function(
    x,
    reference = NULL,
    flow_var = c("NES", "delta_NES", "leading_edge_size"),
    ...
) {
    flow_var <- match.arg(flow_var)

    if (is.list(x) && !inherits(x, c("nseaResult", "mnseaResult"))) {
        if (is.null(names(x))) {
            stop("A list input to mechanismflow() must be named.")
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

    # Stable mechanism-state order.
    state_levels <- c(
        "conserved", "context_specific", "rewired", "contradictory", "unknown"
    )
    state_levels <- intersect(state_levels, unique(df$mechanism_class))
    if (length(state_levels) == 0) {
        state_levels <- "unknown"
    }
    df$mechanism_class <- factor(df$mechanism_class, levels = state_levels)
    df$state_num <- as.integer(df$mechanism_class)

    flow_raw <- as.numeric(df[[flow_var]])
    flow_values <- abs(flow_raw)
    flow_values[is.na(flow_values)] <- 1
    df$flow_size <- flow_values

    p <- ggplot(
        df,
        aes(
            x = .data$context,
            y = .data$state_num,
            group = .data$ID,
            color = .data$Description,
            linewidth = .data$flow_size
        )
    ) +
        geom_line(alpha = 0.8) +
        geom_point(
            aes(size = .data$flow_size, shape = .data$mechanism_class),
            alpha = 0.9
        ) +
        scale_y_continuous(
            breaks = seq_along(state_levels),
            labels = state_levels
        ) +
        scale_size_continuous(name = "Flow magnitude") +
        ggplot2::scale_color_discrete(name = "Pathway") +
        labs(x = "Context", y = "Mechanism state", shape = "State") +
        theme_bw()

    p
}
