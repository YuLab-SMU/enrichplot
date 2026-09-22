#' Plotting utility functions for enrichplot package
#'
#' This file contains plotting and visualization helper functions for enrichplot

#' Automatically split barplot or dotplot into several facets
#'
#' @param by one of 'row' or 'column'
#' @param scales whether 'fixed' or 'free'
#' @param levels set facet levels
#' @return a ggplot object
#' @export
autofacet <- function(by = 'row', scales = "free", levels = NULL) {
    structure(
        list(by = by, scales = scales, levels = levels),
        class = "autofacet"
    )
}

#' Internal plot function for plotting compareClusterResult
#'
#' @param clProf.reshape.df data frame of compareCluster result
#' @param x x variable
#' @param type one of dot and bar
#' @param by one of percentage and count
#' @param title graph title
#' @param font.size graph font size
#' @param colorBy one of pvalue or p.adjust
#' @param width width of the plot
#' @return ggplot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 %+%
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_color_continuous
#' @importFrom ggplot2 guide_colorbar
#' @author Guangchuang Yu <https://yulab-smu.top>
plotting.clusterProfile <- function(
    clProf.reshape.df,
    x = ~Cluster,
    type = "dot",
    colorBy = "p.adjust",
    by = "geneRatio",
    title = "",
    font.size = 12,
    width = NULL
) {
    if (type == "bar") {
        yvar <- switch(by,
            percentage = "Percentage",
            rowPercentage = "Percentage",
            count = "Count",
            geneRatio = "GeneRatio",
            by
        )
        p <- ggplot(
            clProf.reshape.df,
            aes(
                x = !!sym("Description"),
                y = !!sym(yvar),
                fill = !!sym("Cluster")
            )
        )
        if (is.null(width)) {
            p <- p + geom_col()
        } else {
            p <- p + geom_col(width = width)
        }
        p <- p + coord_flip()
    }

    p <- p + xlab("") + ylab("") + ggtitle(title) + theme_dose(font.size)
    return(p)
}

#' Get the distance of the label
#'
#' @param dimension one of 1 and 2
#' @param label_location label_location
#' @return distance matrix
#' @noRd
get_label_diss <- function(dimension, label_location) {
    nn <- nrow(label_location)
    label_dis <- matrix(NA, nrow = nn, ncol = nn)
    colnames(label_dis) <- rownames(label_dis) <- label_location$label

    # Vectorized computation using outer
    vals <- label_location[[dimension]]
    label_dis <- outer(vals, vals, `-`)
    colnames(label_dis) <- rownames(label_dis) <- label_location$label

    # Convert to long format
    label_diss <- reshape2::melt(label_dis)
    label_diss <- label_diss[label_diss[, 1] != label_diss[, 2], ]
    label_diss <- label_diss[!is.na(label_diss[, 3]), ]
    label_diss[, 1] <- as.character(label_diss[, 1])
    label_diss[, 2] <- as.character(label_diss[, 2])
    return(label_diss)
}

#' Default labeller function
#'
#' Default labeling function that uses the
#' internal string wrapping function `yulab.utils::str_wrap`
#' @noRd
#' @importFrom yulab.utils str_wrap
default_labeller <- function(n) {
    fun <- function(str) {
        str <- gsub("_", " ", str)
        yulab.utils::str_wrap(str, n)
    }

    structure(fun, class = "labeller")
}

#' Get segment.size value for ggrepel
#'
#' @param default default value of ggrepel.segment.size
#' @return segment size value
#' @noRd
get_ggrepel_segsize <- function(default = 0.2) {
    getOption("ggrepel.segment.size", default = default)
}


#' Get parameter change message
#'
#' @param parameter parameter name
#' @param params_df parameter data frame
#' @return warning message
#' @noRd
get_param_change_message <- function(parameter, params_df) {
    paste0(
        "Use '",
        params_df[parameter, "listname"],
        " = list(",
        params_df[parameter, "present"],
        " = your_value)' instead of '",
        params_df[parameter, "original"],
        "'"
    )
}
