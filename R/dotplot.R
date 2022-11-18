##' @rdname dotplot
##' @exportMethod dotplot
##' @author guangchuang yu
setMethod("dotplot", signature(object = "enrichResult"),
          function(object, x = "GeneRatio", color = "p.adjust",
                   showCategory=10, size = NULL,
                   split = NULL, font.size=12, title = "",
                   orderBy="x", label_format = 30, ...) {
              dotplot.enrichResult(object = object, x = x, color = color,
                                   showCategory = showCategory,
                                   size = size, split = split,
                                   font.size = font.size,
                                   title = title, orderBy = orderBy,
                                   label_format = label_format, ...)
          })

##' @rdname dotplot
##' @exportMethod dotplot
setMethod("dotplot", signature(object = "gseaResult"),
          function(object, x = "GeneRatio", color = "p.adjust", showCategory=10,
                   size = NULL, split = NULL, font.size=12, title = "",
                   orderBy="x", label_format = 30, ...) {
                if (color == "NES") {
                    NES <- TRUE
                    color <- "p.adjust"
                } else {
                    NES <- FALSE
                }
                p <- dotplot.enrichResult(object = object, x = x, color = color,
                        showCategory = showCategory,
                        size = size, split = split,
                        font.size = font.size,
                        title = title, orderBy = orderBy,
                        label_format = label_format, ...)
                
                if (NES) {
                    p <- suppressMessages(p + aes_(color=~NES) + 
                        scale_color_continuous(low="blue", high="red", name = "NES")
                    )
                }
                return(p)
          })

##' @rdname dotplot
##' @aliases dotplot,compareClusterResult,ANY-method
##' @exportMethod dotplot
setMethod("dotplot", signature(object="compareClusterResult"),
          function(object,
                   x = "Cluster",
                   color ="p.adjust",
                   showCategory=5,
                   split=NULL,
                   font.size=12,
                   title="",
                   by="geneRatio",
                   size=NULL,
                   includeAll=TRUE,
                   label_format = 30,
                   ...
                   ) {
              dotplot.compareClusterResult(object, x=x, colorBy = color,
                                           showCategory = showCategory, by = by,
                                           size = size, includeAll = includeAll,
                                           split = split, font.size = font.size,
                                           title = title, label_format = label_format,
                                           ...)
})

##' @rdname dotplot
##' @exportMethod dotplot
##' @aliases dotplot,enrichResultList,ANY-method
##' @author guangchuang yu
setMethod("dotplot", signature(object = "enrichResultList"),
          function(object, x = "GeneRatio", color = "p.adjust",
                   showCategory=10, size = NULL,
                   split = NULL, font.size=12, title = "",
                   orderBy="x", label_format = 30, ...) {
              dotplot.enrichResult(object = object, x = x, color = color,
                                   showCategory = showCategory,
                                   size = size, split = split,
                                   font.size = font.size,
                                   title = title, orderBy = orderBy,
                                   label_format = label_format, ...)
})

##' @rdname dotplot
##' @exportMethod dotplot
##' @aliases dotplot,gseaResultList,ANY-method
setMethod("dotplot", signature(object = "gseaResultList"),
          function(object, x = "GeneRatio", color = "p.adjust",
                   showCategory=10, size = NULL,
                   split = NULL, font.size=12, title = "",
                   orderBy="x", label_format = 30, ...) {
              if (color == "NES") {
                      NES <- TRUE
                      color <- "p.adjust"
              } else {
                      NES <- FALSE
              }
              p <- dotplot.enrichResult(object = object, x = x, color = color,
                      showCategory = showCategory,
                      size = size, split = split,
                      font.size = font.size,
                      title = title, orderBy = orderBy,
                      label_format = label_format, ...)
              
              if (NES) {
                  p <- suppressMessages(p + aes_(color=~NES) + 
                      scale_color_continuous(low="blue", high="red", name = "NES")
                  )
              }
              return(p)
})

##' @rdname dotplot
##' @param x variable for x-axis, one of 'GeneRatio' and 'Count'
##' @param color variable that used to color enriched terms,
##'              e.g. 'pvalue', 'p.adjust' or 'qvalue'
##' @param showCategory A number or a list of terms. If it is a number,
##' the first n terms will be displayed. If it is a list of terms,
##' the selected terms will be displayed.
##' @param size variable that used to scale the sizes of categories,
##' one of "geneRatio", "Percentage" and "count"
##' @param split separate result by 'category' variable
##' @param font.size font size
##' @param title plot title
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' by default wraps names longer that 30 characters
##' @param orderBy The order of the Y-axis
##' @param decreasing logical. Should the orderBy order be increasing or decreasing?
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
##' @importFrom ggplot2 guides
##' @importFrom ggplot2 guide_legend
##' @importFrom methods is
dotplot.enrichResult <- function(object, x = "geneRatio", color = "p.adjust",
                             showCategory=10, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x",
                             label_format = 30, decreasing=TRUE) {

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
    
    if (inherits(object, c("enrichResultList", "gseaResultList"))) {
        ldf <- lapply(object, fortify, showCategory=showCategory, split=split)
        df <- dplyr::bind_rows(ldf, .id="category")
        df$category <- factor(df$category, levels=names(object))
    } else {
        df <- fortify(object, showCategory = showCategory, split=split)
        ## already parsed in fortify
        ## df$GeneRatio <- parse_ratio(df$GeneRatio)
    }

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
        scale_size(range=c(3, 8)) +
        guides(size  = guide_legend(order = 1),
               color = guide_colorbar(order = 2))
}


##' @rdname dotplot
##' @param object compareClusterResult object
##' @param by one of "geneRatio", "Percentage" and "count"
##' @param split ONTOLOGY or NULL
##' @param includeAll logical
##' @param font.size font size
##' @param title figure title
##' @param group a logical value, whether to connect the
##' nodes of the same group with wires.
##' @param shape a logical value, whether to use nodes of
##' different shapes to distinguish the group it belongs to
##' @param colorBy variable that used to color enriched terms,
##' e.g. 'pvalue', 'p.adjust' or 'qvalue'
##' @importFrom ggplot2 facet_grid
dotplot.compareClusterResult <- function(object, x= "Cluster", colorBy="p.adjust",
                                         showCategory=5, by="geneRatio", size="geneRatio",
                                         split=NULL, includeAll=TRUE,
                                         font.size=12, title="", label_format = 30,
                                         group = FALSE, shape = FALSE) {
    color <- NULL
    if (is.null(size)) size <- by ## by may deprecated in future release

    df <- fortify(object, showCategory=showCategory, by=size,
                  includeAll=includeAll, split=split)

    if (by != "geneRatio")
        df$GeneRatio <- parse_ratio(df$GeneRatio)
    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }
    if (size %in% c("rowPercentage", "count", "geneRatio")) {
        by2 <- switch(size, rowPercentage = "Percentage",
                            count         = "Count",
                            geneRatio     = "GeneRatio")
    } else {
        by2 <- size
    }

    p <- ggplot(df, aes_string(x = x, y = "Description", size = by2)) +
        scale_y_discrete(labels = label_func)

    ## show multiply GO enrichment result in separate panels
    #if ("ONTOLOGY" %in% colnames(df) && length(unique(df$ONTOLOGY)) > 1){
    #    p = p + facet_grid(
    #        ONTOLOGY ~ .,
    #        scales = "free",
    #        space = "free"
    #    )
    #}

    if (group) {
        p <- p + geom_line(aes_string(color = "Cluster", group = "Cluster"), size=.3) +
          ggnewscale::new_scale_colour()
    }

    if (shape) {
        ggstar <- "ggstar"
        require(ggstar, character.only=TRUE)
        # p <- p + ggsymbol::geom_symbol(aes_string(symbol = "Cluster", fill = colorBy)) +
        p <- p + ggstar::geom_star(aes_string(starshape="Cluster", fill=colorBy)) +
            scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))
    }  else {
        p <- p +  geom_point(aes_string(color = colorBy))
    }

    p + scale_color_continuous(low="red", high="blue",
                    guide=guide_colorbar(reverse=TRUE)) +
        ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) +
        scale_size_continuous(range=c(3, 8)) +
        guides(size  = guide_legend(order = 1),
                color = guide_colorbar(order = 2))
}
