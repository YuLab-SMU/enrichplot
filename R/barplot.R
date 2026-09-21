#' Barplot of enrichResult
#'
#' Barplot of enrichResult
#'
#' @importFrom graphics barplot
#' @importFrom ggplot2 %+%
#' @importFrom ggplot2 scale_fill_continuous
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_y_discrete
#' @title barplot
#' @param height enrichResult object
#' @param x one of 'Count' and 'GeneRatio'
#' @param color one of 'pvalue', 'p.adjust' and 'qvalue'
#' @param showCategory number of categories to display or a vector of terms.
#' @param font.size font size
#' @param title plot title
#' @param label_format a numeric value sets wrap length, alternatively a
#' custom function to format axis labels.
#' by default wraps names longer than 30 characters
#' @param ... additional parameters
#' @method barplot enrichResult
#' @export
#' @return ggplot object
#' @examples
#' \dontrun{
#' library(DOSE)
#' data(geneList)
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' barplot(x)
#' # use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
#' barplot(x, showCategory = 10)
#' categories <- c("urinary bladder cancer", "bronchiolitis obliterans",
#'                "aortic aneurysm", "esophageal cancer")
#' barplot(x, showCategory = categories)
#' }
barplot.enrichResult <- function(
    height,
    x = "Count",
    color = "p.adjust",
    showCategory = 8,
    font.size = 12,
    title = "",
    label_format = 30,
    ...
) {
    barplot_internal(
        height = height,
        x = x,
        color = color,
        showCategory = showCategory,
        font.size = font.size,
        title = title,
        label_format = label_format,
        fortify_fun = fortify.enrichResult,
        ...
    )
}

#' @method barplot gseaResult
#' @export
barplot.gseaResult <- function(
    height,
    x = "Count",
    color = "p.adjust",
    showCategory = 8,
    font.size = 12,
    title = "",
    label_format = 30,
    ...
) {
    barplot_internal(
        height = height,
        x = x,
        color = color,
        showCategory = showCategory,
        font.size = font.size,
        title = title,
        label_format = label_format,
        fortify_fun = fortify.gseaResult,
        ...
    )
}

barplot_internal <- function(
    height,
    x = "Count",
    color = "p.adjust",
    showCategory = 8,
    font.size = 12,
    title = "",
    label_format = 30,
    fortify_fun = fortify.enrichResult,
    ...
) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an enrichment result object.
    object <- height

    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    x <- normalize_measure_var(x)

    dots <- list(...)
    supported_params <- c("order", "drop", "split")
    fortify_params <- dots[names(dots) %in% supported_params]
    geom_width <- dots[["width"]]

    # Create the call to fortify without passing ... directly
    # This prevents ggplot2 from checking for unused parameters
    fortify_args <- list(
        model = object,
        showCategory = showCategory,
        by = x
    )

    # Add supported parameters explicitly
    if ("order" %in% names(fortify_params)) {
        fortify_args$order <- fortify_params$order
    }
    if ("drop" %in% names(fortify_params)) {
        fortify_args$drop <- fortify_params$drop
    }
    if ("split" %in% names(fortify_params)) {
        fortify_args$split <- fortify_params$split
    }

    df <- do.call(fortify_fun, fortify_args)

    if (colorBy %in% colnames(df)) {
        p <- ggplot(
            df,
            aes(
                x = .data[[x]],
                y = .data[["Description"]],
                fill = .data[[colorBy]]
            )
        ) +
            theme_dose(font.size) +
            set_enrichplot_color(type = "fill", name = color)
    } else {
        p <- ggplot(
            df,
            aes(
                x = .data[[x]],
                y = .data[["Description"]],
                fill = .data[["Description"]]
            )
        ) +
            theme_dose(font.size) +
            theme(legend.position = "none")
    }

    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    if (is.null(geom_width)) {
        p <- p + geom_col()
    } else {
        p <- p + geom_col(width = geom_width)
    }

    p +
        scale_y_discrete(labels = label_func) +
        ggtitle(title) +
        ylab(NULL) # + xlab(NULL)
}

#' @method barplot compareClusterResult
#' @export
barplot.compareClusterResult <- function(
    height,
    color = "p.adjust",
    showCategory = 5,
    by = "geneRatio",
    includeAll = TRUE,
    font.size = 12,
    title = "",
    ...
) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an compareClusterResult object.
    df <- fortify(
        height,
        showCategory = showCategory,
        by = by,
        includeAll = includeAll
    )
    geom_width <- list(...)[["width"]]
    plotting.clusterProfile(
        df,
        type = "bar",
        colorBy = color,
        by = by,
        title = title,
        font.size = font.size,
        width = geom_width
    )
}
