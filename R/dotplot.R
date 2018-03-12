##' @rdname dotplot
##' @param x variable for x-axis, one of 'geneRatio' or 'Count'
##' @param color variable that used to color enriched terms,
##'              e.g. pvalue, p.adjust or qvalue
##' @param showCategory number of enriched terms to display
##' @param split separate result by 'category' variable
##' @param font.size font size
##' @param title plot title
##' @importClassesFrom DOSE enrichResult
##' @exportMethod dotplot
##' @author guangchuang yu
setMethod("dotplot", signature(object = "enrichResult"),
          function(object, x = "geneRatio", color = "p.adjust", showCategory=10, split = NULL, font.size=12, title = "", ...) {
              dotplot_internal(object, x, color, showCategory, split, font.size, title, ...)
          })

##' @rdname dotplot
##' @importClassesFrom DOSE gseaResult
##' @exportMethod dotplot
setMethod("dotplot", signature(object = "gseaResult"),
          function(object, x = "geneRatio", color = "p.adjust", showCategory=10, split = NULL, font.size=12, title = "", ...) {
              dotplot_internal(object, x, color, showCategory, split, font.size, title, ...)
          })



##' @importFrom ggplot2 fortify
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 scale_color_gradient
##' @importFrom ggplot2 scale_color_continuous
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
dotplot_internal <- function(object, x = "geneRatio", color = "p.adjust", showCategory=10, split = NULL,
                             font.size=12, title = "", orderBy="GeneRatio", decreasing=TRUE) {

    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
        size <- "Count"
    } else if (x == "count" || x == "Count") {
        x <- "Count"
        size <- "GeneRatio"
    } else {
        stop("x should be geneRatio or count...")
    }
    df <- fortify(object, showCategory = showCategory, split=split)
    ## already parsed in fortify
    ## df$GeneRatio <- parse_ratio(df$GeneRatio)

    if (!orderBy %in% colnames(df)) {
        message('wrong orderBy parameter; set to default `orderBy = "GeneRatio"`')
        orderBy <- "GeneRatio"
    }

    idx <- order(df[[orderBy]], decreasing = decreasing)
    df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
        geom_point() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        ylab(NULL) + ggtitle(title) + theme_dose(font.size) + scale_size(range=c(3, 8))

}

