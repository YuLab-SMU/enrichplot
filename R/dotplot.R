#' @rdname dotplot
#' @exportMethod dotplot
#' @author Guangchuang Yu
setMethod(
    "dotplot",
    signature(object = "enrichResult"),
    function(
        object,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        orderBy = "x",
        label_format = 30,
        ...
    ) {
        dotplot.enrichResult(
            object = object,
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            orderBy = orderBy,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname dotplot
#' @exportMethod dotplot
setMethod(
    "dotplot",
    signature(object = "gseaResult"),
    function(
        object,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        orderBy = "x",
        label_format = 30,
        ...
    ) {
        if (color == "NES") {
            NES <- TRUE
            color <- "p.adjust"
        } else {
            NES <- FALSE
        }
        p <- dotplot.enrichResult(
            object = object,
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            orderBy = orderBy,
            label_format = label_format,
            ...
        )

        if (NES) {
            p <- suppressMessages(
                p +
                    aes(fill = .data$NES) +
                    set_enrichplot_color(
                        type = "fill",
                        name = "NES"
                    )
            )
        }
        return(p)
    }
)

#' @rdname dotplot
#' @exportMethod dotplot
setMethod(
    "dotplot",
    signature(object = "mnseaResult"),
    function(
        object,
        x = "share",
        color = "contribution",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        orderBy = "x",
        label_format = 30,
        layer = NULL,
        ...
    ) {
        dotplot.mnseaResult(
            object = object,
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            orderBy = orderBy,
            label_format = label_format,
            layer = layer,
            ...
        )
    }
)

#' @rdname dotplot
#' @aliases dotplot,compareClusterResult,ANY-method
#' @exportMethod dotplot
setMethod(
    "dotplot",
    signature(object = "compareClusterResult"),
    function(
        object,
        x = "Cluster",
        color = "p.adjust",
        showCategory = 5,
        split = NULL,
        font.size = 12,
        title = "",
        by = "geneRatio",
        size = NULL,
        includeAll = TRUE,
        label_format = 30,
        ...
    ) {
        dotplot.compareClusterResult(
            object,
            x = x,
            colorBy = color,
            showCategory = showCategory,
            by = by,
            size = size,
            includeAll = includeAll,
            split = split,
            font.size = font.size,
            title = title,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname dotplot
#' @exportMethod dotplot
#' @aliases dotplot,enrichResultList,ANY-method
#' @author Guangchuang Yu
setMethod(
    "dotplot",
    signature(object = "enrichResultList"),
    function(
        object,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        orderBy = "x",
        label_format = 30,
        ...
    ) {
        dotplot.enrichResult(
            object = object,
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            orderBy = orderBy,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname dotplot
#' @exportMethod dotplot
#' @aliases dotplot,gseaResultList,ANY-method
setMethod(
    "dotplot",
    signature(object = "gseaResultList"),
    function(
        object,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        orderBy = "x",
        label_format = 30,
        ...
    ) {
        if (color == "NES") {
            NES <- TRUE
            color <- "p.adjust"
        } else {
            NES <- FALSE
        }
        p <- dotplot.enrichResult(
            object = object,
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            orderBy = orderBy,
            label_format = label_format,
            ...
        )

        if (NES) {
            p <- suppressMessages(
                p +
                    aes(fill = .data$NES) +
                    set_enrichplot_color(type = "fill", name = "NES")
            )
        }
        return(p)
    }
)

#' @rdname dotplot
#' @param x variable for x-axis, one of 'GeneRatio' and 'Count'
#' @param color variable used to color enriched terms,
#'              e.g. 'pvalue', 'p.adjust' or 'qvalue'
#' @param showCategory number of categories to display or a vector of terms.
#' @param size variable used to scale the sizes of categories,
#' one of "geneRatio", "Percentage" and "count"
#' @param split separate result by 'category' variable
#' @param font.size font size
#' @param title plot title
#' @param label_format a numeric value sets wrap length, alternatively a
#' custom function to format axis labels.
#' by default wraps names longer than 30 characters
#' @param orderBy The order of the Y-axis
#' @param decreasing logical. Should the orderBy order be increasing or decreasing?
#' @param layer optional layer or layers to retain for `mnseaResult`.
#' @importFrom ggplot2 fortify
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 scale_color_continuous
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom methods is
dotplot.enrichResult <- function(
    object,
    x = "geneRatio",
    color = "p.adjust",
    showCategory = 10,
    size = NULL,
    split = NULL,
    font.size = 12,
    title = "",
    orderBy = "x",
    label_format = 30,
    decreasing = TRUE
) {
    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    .formula_expr <- NULL
    x <- normalize_measure_var(x)
    if (x == "GeneRatio") {
        if (is.null(size)) {
            size <- "Count"
        }
    } else if (x == "Count") {
        if (is.null(size)) {
            size <- "GeneRatio"
        }
    } else if (is(x, "formula")) {
        .formula_expr <- rlang::f_rhs(x)
        x <- "x"
        if (is.null(size)) {
            size <- "Count"
        }
    } else {
        ## message("invalid x, setting to 'GeneRatio' by default")
        ## x <- "GeneRatio"
        ## size <- "Count"
        if (is.null(size)) {
            size <- "Count"
        }
    }

    if (inherits(object, c("enrichResultList", "gseaResultList"))) {
        ldf <- lapply(
            object,
            fortify,
            showCategory = showCategory,
            split = split
        )
        df <- dplyr::bind_rows(ldf, .id = "category")
        df$category <- factor(df$category, levels = names(object))
    } else {
        df <- fortify(object, showCategory = showCategory, split = split)
        ## already parsed in fortify
        ## df$GeneRatio <- parse_ratio(df$GeneRatio)
    }

    if (!is.null(.formula_expr)) {
        df$x <- rlang::eval_tidy(.formula_expr, data = df)
    }

    if (orderBy != 'x' && !orderBy %in% colnames(df)) {
        message('wrong orderBy parameter; set to default `orderBy = "x"`')
        orderBy <- "x"
    }

    if (orderBy == "x") {
        df <- dplyr::mutate(df, x = .data[[x]])
    }

    label_func <- .label_format(label_format)

    idx <- order(df[[orderBy]], decreasing = decreasing)

    df$Description <- factor(
        df$Description,
        levels = rev(unique(df$Description[idx]))
    )

    # Use internal helper function for common plotting logic
    p <- .dotplot_internal(
        df = df,
        x = x,
        size = size,
        colorBy = colorBy,
        color = color,
        label_func = label_func,
        font.size = font.size,
        title = title
    )
    
    return(p)
}

.label_format <- function(label_format = 30) {
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    return(label_func)
}

dotplot.mnseaResult <- function(
    object,
    x = "share",
    color = "contribution",
    showCategory = 10,
    size = NULL,
    split = NULL,
    font.size = 12,
    title = "",
    orderBy = "x",
    label_format = 30,
    layer = NULL,
    decreasing = TRUE
) {
    x <- match.arg(x, c("share", "contribution"))
    colorBy <- match.arg(
        color,
        c("contribution", "NES", "pvalue", "p.adjust", "qvalue")
    )
    if (is.null(size)) {
        size <- "n_feature"
    }

    fortify_by <- if (x == "contribution") {
        "contribution"
    } else if (colorBy %in% c("contribution", "NES")) {
        colorBy
    } else {
        "p.adjust"
    }

    df <- fortify(
        object,
        showCategory = showCategory,
        by = fortify_by,
        level = "pathway",
        layer = layer
    )

    if (nrow(df) == 0) {
        yulab.utils::yulab_abort("No mnsea contribution data available for plotting.")
    }

    if (orderBy != "x" && !orderBy %in% colnames(df)) {
        message('wrong orderBy parameter; set to default `orderBy = "x"`')
        orderBy <- "x"
    }

    if (orderBy == "x") {
        df <- dplyr::mutate(df, x = .data[[x]])
    }

    label_func <- .label_format(label_format)
    idx <- order(df[[orderBy]], decreasing = decreasing)
    df$Description <- factor(
        as.character(df$Description),
        levels = rev(unique(as.character(df$Description[idx])))
    )
    df <- filter_mnsea_layers(df, layer = layer)
    if (nrow(df) == 0) {
        yulab.utils::yulab_abort("No mnsea contribution data available for the selected `layer`.")
    }

    color_scale <- switch(
        colorBy,
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
        list(
            colors = get_enrichplot_color(2),
            transform = "log10",
            reverse = TRUE
        )
    )

    p <- .dotplot_internal(
        df = df,
        x = x,
        size = size,
        colorBy = colorBy,
        color = mnsea_plot_label(colorBy),
        label_func = label_func,
        font.size = font.size,
        title = title,
        size_name = mnsea_plot_label(size),
        color_colors = color_scale$colors,
        color_transform = color_scale$transform,
        color_reverse = color_scale$reverse
    )

    if (!is.null(split) && split %in% colnames(df)) {
        p <- p +
            facet_grid(
                .data[[split]] ~ .,
                scales = "free_y",
                space = "free_y",
                switch = "y"
            )
    }

    p
}

#' Internal helper function for dotplot construction
#' @param df data frame containing the plot data
#' @param x x-axis variable name
#' @param size size variable name
#' @param colorBy color variable name
#' @param color color parameter name for legend
#' @param label_func function for formatting labels
#' @param font.size font size for theme
#' @param title plot title
#' @param size_range range for size scaling, default c(3, 8)
#' @param shape_point whether to use enrichplot_point_shape, default TRUE
#' @param color_colors color palette used for fill mapping
#' @param color_transform transform applied to the fill scale
#' @param color_reverse whether to reverse the color legend
#' @return ggplot object with enrichplotDot class
#' @noRd
.dotplot_internal <- function(
    df,
    x,
    size,
    colorBy,
    color,
    label_func,
    font.size,
    title,
    size_range = c(3, 8),
    size_name = ggplot2::waiver(),
    shape_point = TRUE,
    color_colors = get_enrichplot_color(2),
    color_transform = "log10",
    color_reverse = TRUE
) {
    p <- ggplot(
        df,
        aes(
            x = .data[[x]],
            y = .data[["Description"]],
            size = .data[[size]],
            fill = .data[[colorBy]]
        )
    ) +
        geom_point() +
        set_enrichplot_color(
            colors = color_colors,
            type = "fill",
            name = color,
            transform = color_transform,
            reverse = color_reverse
        ) +
        scale_y_discrete(labels = label_func) +
        ylab(NULL) +
        ggtitle(title) +
        theme_dose(font.size)
    
    if (shape_point) {
        p <- p + aes(shape = I(enrichplot_point_shape))
    }
    
    # Apply size scaling
    if (size == "Count") {
        # For Count, use pretty breaks
        size_break <- pretty(df[[size]], n = 4)
        p <- p + scale_size(range = size_range, breaks = size_break, name = size_name)
    } else {
        p <- p + scale_size(range = size_range, name = size_name)
    }
    
    class(p) <- c("enrichplotDot", class(p))
    return(p)
}


#' @rdname dotplot
#' @param object compareClusterResult object
#' @param by one of "geneRatio", "Percentage" and "count"
#' @param split apply `showCategory` to each category specified by the 'split', e.g., "ONTOLOGY", "category" and "intersect".  Default is NULL and do nothing
#' @param includeAll logical value
#' @param font.size font size
#' @param title figure title
#' @param group a logical value, whether to connect the
#' nodes of the same group with wires.
#' @param shape a logical value, whether to use nodes of
#' different shapes to distinguish the group it belongs to
#' @param facet apply `facet_grid` to the plot by specified variable, e.g., "ONTOLOGY", "category" and "intersect".
#' @param strip_width width of strip text (facet label).
#' @param colorBy variable used to color enriched terms,
#' e.g. 'pvalue', 'p.adjust' or 'qvalue'
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom rlang check_installed
dotplot.compareClusterResult <- function(
    object,
    x = "Cluster",
    colorBy = "p.adjust",
    showCategory = 5,
    by = "geneRatio",
    size = "geneRatio",
    split = NULL,
    includeAll = TRUE,
    font.size = 12,
    title = "",
    label_format = 30,
    group = FALSE,
    shape = FALSE,
    facet = NULL,
    strip_width = 15
) {
    color <- NULL
    if (is.null(size)) {
        size <- by
    } ## by may deprecated in future release

    if (!is.null(facet) && facet == "intersect") {
        object <- append_intersect(object)
    }

    df <- fortify(
        object,
        showCategory = showCategory,
        by = size,
        includeAll = includeAll,
        split = split
    )

    # if (by != "geneRatio")
    #    df$GeneRatio <- parse_ratio(df$GeneRatio)
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    by2 <- normalize_measure_var(size, include_percentage = TRUE)

    # Use internal helper function for base plot, but without shape_point for flexibility
    p <- .dotplot_internal(
        df = df,
        x = x,
        size = by2,
        colorBy = colorBy,
        color = colorBy,
        label_func = label_func,
        font.size = font.size,
        title = title,
        shape_point = FALSE  # We'll handle shape manually
    )
    
    # Add group connections if requested
    if (group) {
        p <- p +
            geom_line(
                aes(
                    x = .data[[x]],
                    y = .data[["Description"]],
                    color = .data$Cluster,
                    group = .data$Cluster
                ),
                linewidth = .3,
                inherit.aes = FALSE
            ) +
            ggnewscale::new_scale_colour()
    }
    
    # Handle shape variations
    if (shape) {
        require_suggested('ggstar', 'for `dotplot()` with `shape = TRUE`.')
        ggstar <- "ggstar"
        require(ggstar, character.only = TRUE)
        # Replace the base geom_point with ggstar
        p$layers <- p$layers[-1]  # Remove the first geom_point layer
        p <- p +
            ggstar::geom_star(aes(
                starshape = .data$Cluster,
                fill = .data[[colorBy]]
            )) +
            set_enrichplot_color(type = "fill", transform = 'log10')
    } else {
        # Add standard point with shape
        p <- p + aes(shape = I(enrichplot_point_shape))
    }
    
    # Add facet if requested
    if (!is.null(facet)) {
        p <- p +
            facet_grid(
                .data[[facet]] ~ .,
                scales = "free",
                space = 'free',
                switch = 'y',
                labeller = ggplot2::label_wrap_gen(strip_width)
            ) +
            theme(strip.text = element_text(size = 14))
    }
    
    return(p)
}


append_intersect <- function(x) {
    if (!inherits(x, 'compareClusterResult')) {
        stop("x should be a compareClusterResult object")
    }

    d <- as.data.frame(x)

    so <- vapply(
        split(d$Cluster, d$Description),
        FUN = paste,
        FUN.VALUE = character(1),
        collapse = " & "
    )

    set_info <- data.frame(
        intersect = so,
        Description = names(so)
    )

    d2 <- merge(d, set_info, by = "Description")
    n <- levels(d2$Cluster)
    cc <- yulab.utils::combinations(length(n))
    lv <- vapply(cc, function(i) paste(n[i], collapse = " & "), character(1))
    d2$intersect <- factor(d2$intersect, levels = lv)

    x@compareClusterResult <- d2

    return(x)
}

#' compare two clusters in the compareClusterResult object
#'
#'
#' @title dotplot2
#' @param object a compareClusterResult object
#' @param x selected variable to visualize in x-axis
#' @param vars selected Clusters to be compared, only length of two is supported
#' @param label to label the Clusters in the plot, should be a named vector
#' @param ... additional parameters passed to dotplot
#' @return a ggplot object
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_blank
#' @importFrom grid arrow
#' @importFrom grid unit
#' @importFrom ggplot2 geom_text
#' @export
#' @author Guangchuang Yu
dotplot2 <- function(
    object,
    x = "FoldEnrichment",
    vars = NULL,
    label = "auto",
    ...
) {
    if (!is(object, 'compareClusterResult')) {
        stop('only compareClusterResult object is supported')
    }

    if (length(vars) != 2) {
        stop("vars should be of length 2.")
    }
    object <- dplyr::filter(object, .data$Cluster %in% vars)
    d <- object@compareClusterResult
    x <- normalize_measure_var(x)
    if (!x %in% colnames(d)) {
        if (identical(x, "FoldEnrichment")) {
            d$FoldEnrichment <- compute_fold_enrichment(d)
        } else {
            yulab.utils::yulab_abort(
                paste0(
                    "`x` column `",
                    x,
                    "` is not available in `compareClusterResult`."
                )
            )
        }
    }
    d[[x]] <- d[[x]] * ifelse(d$Cluster == vars[1], -1, 1)
    object@compareClusterResult <- d
    p <- dotplot(object, x = x, ...)
    p <- p +
        geom_segment(
            aes(
                x = .data[[x]],
                y = .data[["Description"]],
                xend = 0,
                yend = .data$Description
            ),
            linewidth = 1,
            color = 'grey50',
            inherit.aes = FALSE
        ) +
        geom_vline(xintercept = 0, lty = 'dashed') +
        scale_x_continuous(labels = abs)

    p$layers <- rev(p$layers)

    if (is.null(label)) {
        return(p)
    }

    d <- dplyr::group_by(object@compareClusterResult, .data$Cluster) |>
        dplyr::summarise(mid = max(abs(.data[[x]])) * sign(max(.data[[x]])) / 2)

    if (!identical(label, "auto")) {
        if (is.null(names(label))) {
            yulab.utils::yulab_abort(
                "`label` should be a named vector keyed by cluster names."
            )
        }
        d$Cluster <- unname(label[as.character(d$Cluster)])
    }

    p +
        geom_text(
            aes(x = .data$mid, y = .4, label = .data$Cluster),
            data = d,
            inherit.aes = FALSE,
            size = 5
        ) +
        geom_segment(
            aes(x = 0, xend = .data$mid, y = .2, yend = .2),
            data = d,
            inherit.aes = FALSE,
            arrow = arrow(length = unit(0.30, "cm"), type = "closed")
        ) +
        geom_blank(
            data = data.frame(y = 0),
            aes(y = .data$y),
            inherit.aes = FALSE
        )
}
