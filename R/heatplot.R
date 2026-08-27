#' @rdname heatplot
#' @exportMethod heatplot
setMethod(
    "heatplot",
    signature(x = "enrichResult"),
    function(
        x,
        showCategory = 30,
        showTop = NULL,
        symbol = "rect",
        foldChange = NULL,
        pvalue = NULL,
        label_format = 30,
        pathway_id = NULL,
        layer = NULL,
        value = c("score", "abs_score", "share", "contribution"),
        ...
    ) {
        heatplot.enrichResult(
            x,
            showCategory = showCategory,
            showTop = showTop,
            symbol = symbol,
            foldChange = foldChange,
            pvalue = pvalue,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname heatplot
#' @exportMethod heatplot
setMethod(
    "heatplot",
    signature(x = "gseaResult"),
    function(
        x,
        showCategory = 30,
        showTop = NULL,
        symbol = "rect",
        foldChange = NULL,
        pvalue = NULL,
        label_format = 30,
        pathway_id = NULL,
        layer = NULL,
        value = c("score", "abs_score", "share", "contribution"),
        ...
    ) {
        heatplot.enrichResult(
            x,
            showCategory = showCategory,
            showTop = showTop,
            symbol = symbol,
            foldChange = foldChange,
            pvalue = pvalue,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname heatplot
#' @exportMethod heatplot
setMethod(
    "heatplot",
    signature(x = "mnseaResult"),
    function(
        x,
        showCategory = 30,
        showTop = NULL,
        symbol = "rect",
        foldChange = NULL,
        pvalue = NULL,
        label_format = 30,
        pathway_id = NULL,
        layer = NULL,
        value = c("score", "abs_score", "share", "contribution"),
        ...
    ) {
        heatplot.mnseaResult(
            x,
            showCategory = showCategory,
            showTop = showTop,
            pathway_id = pathway_id,
            layer = layer,
            value = value,
            label_format = label_format,
            ...
        )
    }
)


#' @rdname heatplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom rlang check_installed
#' @param showTop number of top genes ranked by `abs(foldChange) * frequency`
#' to be shown in the heatmap, default NULL means all genes are shown
#' @param label_format a numeric value sets wrap length, alternatively a
#' custom function to format axis labels.
#' by default wraps names longer than 30 characters
#' @param symbol symbol of the nodes, one of "rect" (the default) or "dot"
#' @param pvalue pvalue of genes
#' @param pathway_id optional pathway ID for pathway-specific `mnseaResult`
#' heatmaps.
#' @param layer optional layer or layers to retain for `mnseaResult` plots.
#' @param value fill value for `mnseaResult`; use `"share"` or
#' `"contribution"` for term-layer heatmaps and `"score"` or `"abs_score"`
#' for pathway-specific feature heatmaps.
#' @author Guangchuang Yu
prepare_heatplot_data <- function(x, showCategory, showTop, foldChange, pvalue) {
    selected <- select_terms(x, showCategory)
    geneSets <- selected$geneSets

    if (!is.null(showTop) && showTop > 0) {
        if (is.null(foldChange)) {
            yulab.utils::yulab_abort(
                "`showTop` requires `foldChange`.",
                class = "missing_foldchange_error"
            )
        }
        nfreq <- table(unlist(geneSets))
        nfc <- nfreq * abs(foldChange[names(nfreq)])
        topgenes <- head(names(sort(nfc, decreasing = TRUE)), showTop)
        geneSets <- lapply(geneSets, function(s) intersect(s, topgenes))
        geneSets <- set_geneSet_labels(geneSets, get_geneSet_labels(selected$geneSets))
    }

    foldChange <- fc_readable(x, foldChange)
    pvalue <- fc_readable(x, pvalue)
    d <- list2df(geneSets)
    d$categoryID <- get_geneSet_labels(geneSets)[as.character(d$categoryID)]
    if (!is.null(foldChange)) {
        d$foldChange <- foldChange[as.character(d$Gene)]
    }

    if (!is.null(pvalue)) {
        d$pvalue <- pvalue[as.character(d$Gene)]
    }

    d
}

prepare_heatplot_mnsea_data <- function(
    x,
    showCategory,
    pathway_id,
    layer,
    showTop,
    value
) {
    if (is.null(pathway_id) && value %in% c("share", "contribution")) {
        df <- fortify(
            x,
            showCategory = showCategory,
            by = value,
            level = "pathway",
            layer = layer
        )
        if (nrow(df) == 0) {
            yulab.utils::yulab_abort("No mnsea pathway contribution data available for plotting.")
        }

        df$axis_label <- factor(
            as.character(df$Description),
            levels = levels(df$Description)
        )
        df$fill_value <- df[[value]]

        return(list(
            data = df,
            y_var = "axis_label",
            fill_name = mnsea_plot_label(value),
            colors = get_enrichplot_color(2),
            reverse = FALSE
        ))
    }

    if (!value %in% c("score", "abs_score")) {
        yulab.utils::yulab_abort(
            "When `pathway_id` is provided, `value` must be `score` or `abs_score`."
        )
    }

    pathway_id <- resolve_mnsea_pathway_id(x, pathway_id, level = "feature")
    df <- fortify_mnsea_contribution(
        x,
        pathway_id = pathway_id,
        level = "feature",
        layer = layer
    )
    if (nrow(df) == 0) {
        yulab.utils::yulab_abort(
            "No mnsea feature contribution data available for the selected `pathway_id`."
        )
    }

    feature_rank <- stats::aggregate(
        df$abs_score,
        by = list(Feature = df$Feature),
        FUN = max,
        na.rm = TRUE
    )
    feature_rank <- feature_rank[order(feature_rank$x, decreasing = TRUE), , drop = FALSE]

    if (!is.null(showTop) && showTop > 0) {
        keep_features <- head(feature_rank$Feature, showTop)
        df <- df[df$Feature %in% keep_features, , drop = FALSE]
        feature_rank <- feature_rank[feature_rank$Feature %in% keep_features, , drop = FALSE]
    }

    df$axis_label <- factor(
        as.character(df$Feature),
        levels = rev(feature_rank$Feature)
    )
    df$fill_value <- df[[value]]

    list(
        data = df,
        y_var = "axis_label",
        fill_name = mnsea_plot_label(value),
        colors = if (value == "score") get_enrichplot_color(3) else get_enrichplot_color(2),
        reverse = FALSE
    )
}

heatplot.enrichResult <- function(
    x,
    showCategory = 30,
    showTop = NULL,
    symbol = "rect",
    foldChange = NULL,
    pvalue = NULL,
    label_format = 30
) {
    symbol <- match.arg(symbol, c("rect", "dot"))
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    d <- prepare_heatplot_data(
        x = x,
        showCategory = showCategory,
        showTop = showTop,
        foldChange = foldChange,
        pvalue = pvalue
    )

    p <- ggplot(d, aes(x = .data$Gene, y = .data$categoryID))

    if (symbol == "rect") {
        p <- p + geom_tile(color = 'white')
    }

    get_dotp <- function(p, foldChange, pvalue) {
        if (is.null(foldChange) && is.null(pvalue)) {
            p <- p +
                geom_point(
                    color = 'black',
                    shape = 21,
                    fill = "black",
                    size = 5
                )
            return(p)
        }
        if (!is.null(foldChange) && !is.null(pvalue)) {
            p <- p + geom_point(color = 'black', shape = 21)
            return(p)
        }

        if (is.null(foldChange)) {
            p <- p + geom_point(color = 'black', shape = 21, fill = "black")
        } else {
            p <- p + geom_point(color = 'black', shape = 21, size = 5)
        }

        return(p)
    }

    # copy from https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
    reverselog_trans <- function(base = exp(1)) {
        trans <- function(x) -log(x, base)

        require_suggested('scales', 'for `heatplot()`.')

        inv <- function(x) base^(-x)
        scales::trans_new(
            paste0("reverselog-", format(base)),
            trans,
            inv,
            scales::log_breaks(base = base),
            domain = c(1e-100, Inf)
        )
    }

    if (symbol == "dot") {
        p <- get_dotp(p, foldChange, pvalue)
        ## only dot need size(pvalue) parameter
        if (!is.null(pvalue)) {
            p <- p +
                aes(size = .data$pvalue) +
                scale_size_continuous(
                    range = c(3, 8),
                    trans = reverselog_trans(10)
                )
        }
    }

    if (!is.null(foldChange)) {
        p <- p +
            aes(fill = !!sym('foldChange')) +
            set_enrichplot_color(
                colors = get_enrichplot_color(3),
                type = "fill",
                reverse = FALSE,
                transform = 'identity'
            )
    }

    p +
        xlab(NULL) +
        ylab(NULL) +
        theme_minimal() +
        scale_y_discrete(labels = label_func) +
        theme(
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 1)
        )
}

heatplot.mnseaResult <- function(
    x,
    showCategory = 10,
    pathway_id = NULL,
    layer = NULL,
    showTop = NULL,
    value = c("score", "abs_score", "share", "contribution"),
    label_format = 30
) {
    value <- match.arg(value)
    label_func <- .label_format(label_format)
    plot_data <- prepare_heatplot_mnsea_data(
        x = x,
        showCategory = showCategory,
        pathway_id = pathway_id,
        layer = layer,
        showTop = showTop,
        value = value
    )
    d <- plot_data$data

    p <- ggplot(
        d,
        aes(
            x = .data$layer,
            y = .data[[plot_data$y_var]],
            fill = .data$fill_value
        )
    ) +
        geom_tile(color = "white") +
        xlab(NULL) +
        ylab(NULL) +
        theme_minimal() +
        scale_y_discrete(labels = label_func) +
        theme(
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    if (value == "score") {
        return(
            p + set_enrichplot_color(
                colors = plot_data$colors,
                type = "fill",
                name = plot_data$fill_name,
                transform = "identity",
                reverse = plot_data$reverse
            )
        )
    }

    p + set_enrichplot_color(
        colors = plot_data$colors,
        type = "fill",
        name = plot_data$fill_name,
        transform = "identity",
        reverse = plot_data$reverse
    )
}
