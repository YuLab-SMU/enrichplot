
##' barplot of enrichResult
##'
##'
##' @importFrom graphics barplot
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 scale_fill_continuous
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
## @S3method barplot enrichResult
##' @title barplot
##' @param height enrichResult object
##' @param x one of 'Count' and 'GeneRatio'
##' @param color one of 'pvalue', 'p.adjust', 'qvalue'
##' @param showCategory number of categories to show
##' @param font.size font size
##' @param title plot title
##' @param ... other parameter, ignored
##' @method barplot enrichResult
##' @export
barplot.enrichResult <- function(height, x="Count", color='p.adjust', showCategory=8, font.size=12, title="", ...) {
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
        p <- ggplot(df, aes_string(x = "Description", y = x, fill = colorBy)) +
            theme_dose(font.size) +
            scale_fill_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE))
    } else {
        p <- ggplot(df, aes_string(x = "Description", y = x, fill = "Description")) +
            theme_dose(font.size) +
            theme(legend.position="none")
    }
    p + geom_bar(stat = "identity") + coord_flip() +
        ggtitle(title) + xlab(NULL) + ylab(NULL)
}
