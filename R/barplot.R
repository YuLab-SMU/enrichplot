
##' barplot of enrichResult
##'
##'
##' @importFrom graphics barplot
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 scale_fill_continuous
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_col
##  @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 scale_y_discrete
##' @title barplot
##' @param height enrichResult object
##' @param x one of 'Count' and 'GeneRatio'
##' @param color one of 'pvalue', 'p.adjust', 'qvalue'
##' @param showCategory number of categories to show
##' @param font.size font size
##' @param title plot title
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' by default wraps names longer that 30 characters
##' @param ... other parameter, ignored
##' @method barplot enrichResult
##' @export
##' @return ggplot object
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' barplot(x)
barplot.enrichResult <- function(height, x="Count", color='p.adjust',
                                 showCategory=8, font.size=12, title="",
                                 label_format=30, ...) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an enrichResult object.
    object <- height

    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
    }
    else if (x == "count" || x == "Count") {
        x <- "Count"
    }

    df <- fortify(object, showCategory=showCategory, by=x, ...)

    if(colorBy %in% colnames(df)) {
        p <- ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
            theme_dose(font.size) +
            scale_fill_continuous(low="red", high="blue", name = color,
                                  guide=guide_colorbar(reverse=TRUE))
    } else {
        p <- ggplot(df, aes_string(x = x, y = "Description",
                                   fill = "Description")) +
            theme_dose(font.size) +
            theme(legend.position="none")
    }

    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }

    p + geom_col() + # geom_bar(stat = "identity") + coord_flip() +
        scale_y_discrete(labels = label_func) +
        ggtitle(title) + xlab(NULL) + ylab(NULL)
}


barplot.compareClusterResult <- function(height, color="p.adjust",
                                         showCategory=5, by="geneRatio",
                                         includeAll=TRUE, font.size=12,
                                         title="", ...) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an compareClusterResult object.
    df <- fortify(height, showCategory=showCategory, by=by,
                  includeAll=includeAll)
    plotting.clusterProfile(df, type="bar", colorBy=color, by=by, title=title,
                            font.size=font.size)
}
