##' @rdname dotplot
##' @param x variable for x-axis, one of 'GeneRatio' or 'Count'
##' @param color variable that used to color enriched terms,
##'              e.g. pvalue, p.adjust or qvalue
##' @param showCategory number of enriched terms to display
##' @param size variable that used to scale the sizes of categories
##' @param split separate result by 'category' variable
##' @param font.size font size
##' @param title plot title
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' by default wraps names longer that 30 characters
##' @importClassesFrom DOSE enrichResult
##' @exportMethod dotplot
##' @author guangchuang yu
setMethod("dotplot", signature(object = "enrichResult"),
          function(object, x = "GeneRatio", color = "p.adjust",
                   showCategory=10, size = NULL,
                   split = NULL, font.size=12, title = "",
                   label_format = 30, ...) {
              dotplot_internal(object, x, color, showCategory, size,
                               split, font.size, title, label_format, ...)
          })

##' @rdname dotplot
##' @importClassesFrom DOSE gseaResult
##' @exportMethod dotplot
setMethod("dotplot", signature(object = "gseaResult"),
          function(object, x = "GeneRatio", color = "p.adjust", showCategory=10,
                   size = NULL, split = NULL, font.size=12, title = "",
                   label_format = 30, ...) {
              dotplot_internal(object, x, color, showCategory, size, split,
                               font.size, title, label_format, ...)
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
##' @importFrom ggplot2 scale_y_discrete
dotplot_internal <- function(object, x = "geneRatio", color = "p.adjust",
                             showCategory=10, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x",
                             decreasing=TRUE, label_format = 30) {

    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
        if (is.null(size))
            size <- "Count"
    } else if (x == "count" || x == "Count") {
        x <- "Count"
        if (is.null(size))
            size <- "GeneRatio"
    } else if (is(x, "formula")) {
        x <- as.character(x)[2]
        if (is.null(size))
            size <- "Count"
    } else {
        ## message("invalid x, setting to 'GeneRatio' by default")
        ## x <- "GeneRatio"
        ## size <- "Count"
        if (is.null(size))
            size  <- "Count"
    }

    df <- fortify(object, showCategory = showCategory, split=split)
    ## already parsed in fortify
    ## df$GeneRatio <- parse_ratio(df$GeneRatio)

    if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
        message('wrong orderBy parameter; set to default `orderBy = "x"`')
        orderBy <- "x"
    }

    if (orderBy == "x") {
        df <- dplyr::mutate(df, x = eval(parse(text=x)))
    }

    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }

    idx <- order(df[[orderBy]], decreasing = decreasing)
    df$Description <- factor(df$Description,
                          levels=rev(unique(df$Description[idx])))
    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
        geom_point() +
        scale_color_continuous(low="red", high="blue", name = color,
            guide=guide_colorbar(reverse=TRUE)) +
        scale_y_discrete(labels = label_func) +
        ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
        scale_size(range=c(3, 8))

}


##' dot plot method
##'
##'
##' @docType methods
##' @title dotplot
##' @rdname dotplot-methods
##' @aliases dotplot,compareClusterResult,ANY-method
##' @param object compareClusterResult object
##' @param x x variable
##' @param color one of pvalue or p.adjust
##' @param showCategory category numbers
##' @param by one of geneRatio, Percentage or count
##' @param split ONTOLOGY or NULL
##' @param includeAll logical
##' @param font.size font size
##' @param title figure title
##' @exportMethod dotplot
setMethod("dotplot", signature(object="compareClusterResult"),
          function(object,
                   x = ~Cluster,
                   color ="p.adjust",
                   showCategory=5,
                   split=NULL,
                   font.size=12,
                   title="",
                   by="geneRatio",
                   includeAll=TRUE
                   ) {
              dotplot.compareClusterResult(object, x=x, colorBy = color,
                                           showCategory = showCategory, by = by,
                                           includeAll = includeAll,
                                           split=split, font.size = font.size,
                                           title = title)
          })


dotplot.compareClusterResult <- function(object, x=~Cluster, colorBy="p.adjust",
                                         showCategory=5, by="geneRatio",
                                         split=NULL, includeAll=TRUE,
                                         font.size=12, title="") {

    df <- fortify(object, showCategory=showCategory, by=by,
                  includeAll=includeAll, split=split)
    plotting.clusterProfile(df, x=x, type="dot", colorBy=colorBy,
                            by=by, title=title, font.size=font.size)
}
