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
                    p <- suppressMessages(p + aes_(fill=~NES) + 
                        # scale_fill_continuous(name = "NES") +
                        set_enrichplot_color(type = "fill", name = "NES")
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
                  p <- suppressMessages(p + aes_(fill=~NES) + 
                      # scale_fill_continuous(name = "NES") + 
                      set_enrichplot_color(type = "fill", name = "NES")
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
    p <- ggplot(df, aes_string(x=x, y="Description", size=size, fill=colorBy)) +
        geom_point() +
        aes(shape = I(enrichplot_point_shape)) + 
        # scale_fill_continuous(name = color) +
        set_enrichplot_color(type = "fill", name = color) + 
        scale_y_discrete(labels = label_func) +
        ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
        scale_size(range=c(3, 8))
    
    class(p) <- c("enrichplotDot", class(p))
    return(p)
}


##' @rdname dotplot
##' @param object compareClusterResult object
##' @param by one of "geneRatio", "Percentage" and "count"
##' @param split apply `showCategory` to each category specified by the 'split', e.g., "ONTOLOGY", "category" and "intersect".  Default is NULL and do nothing
##' @param includeAll logical
##' @param font.size font size
##' @param title figure title
##' @param group a logical value, whether to connect the
##' nodes of the same group with wires.
##' @param shape a logical value, whether to use nodes of
##' different shapes to distinguish the group it belongs to
##' @param facet apply `facet_grid` to the plot by specified variable, e.g., "ONTOLOGY", "category" and "intersect".
##' @param strip_width width of strip text, a.k.a facet label.
##' @param colorBy variable that used to color enriched terms,
##' e.g. 'pvalue', 'p.adjust' or 'qvalue'
##' @importFrom ggplot2 facet_grid
##' @importFrom rlang check_installed  
dotplot.compareClusterResult <- function(object, x= "Cluster", colorBy="p.adjust",
                                         showCategory=5, by="geneRatio", size="geneRatio",
                                         split=NULL, includeAll=TRUE,
                                         font.size=12, title="", label_format = 30,
                                         group = FALSE, shape = FALSE, facet=NULL, strip_width=15) {                                   
    color <- NULL
    if (is.null(size)) size <- by ## by may deprecated in future release

    if (!is.null(facet) && facet == "intersect") {
        object <- append_intersect(object)
    }

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

    if (group) {
        p <- p + geom_line(aes_string(color = "Cluster", group = "Cluster"), size=.3) +
          ggnewscale::new_scale_colour()
    }


    if (shape) {
        check_installed('ggstar', 'for `dotplot()` with `shape = TRUE`.')
        ggstar <- "ggstar"
        require(ggstar, character.only=TRUE)
        # p <- p + ggsymbol::geom_symbol(aes_string(symbol = "Cluster", fill = colorBy)) +
        p <- p + ggstar::geom_star(aes_string(starshape="Cluster", fill=colorBy)) +
            set_enrichplot_color(type = "fill")
    }  else {
        p <- p +  geom_point(aes_string(fill = colorBy)) + 
            aes(shape = I(enrichplot_point_shape))
    }

    p <- p + set_enrichplot_color(type = "fill") +
        ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) +
        scale_size_continuous(range=c(3, 8))


    if (!is.null(facet)) {
        p <- p + facet_grid(.data[[facet]] ~ ., 
                scales = "free", space = 'free',
                switch = 'y',
                labeller = ggplot2::label_wrap_gen(strip_width)) +
            theme(strip.text = element_text(size = 14))
    }

    class(p) <- c("enrichplotDot", class(p))
    
    return(p)
}


append_intersect <- function(x) {
    if (!inherits(x, 'compareClusterResult')) stop("x should be a compareClusterResult object")

    d <- as.data.frame(x)
    # sets <- split(d$Description, d$Cluster)
    # yulab.utils::check_pkg('aplotExtra', 'for upsetplot of compareCluster Result')
    # tidy_main_subsets <- yulab.utils::get_fun_from_pkg('aplotExtra', 'tidy_main_subsets')

    # df <- tidy_main_subsets(sets, order.intersect.by = 'name', nintersects = 10, order.set.by = 'name')
    # nn <- names(sets)
    # df2 <- tidyr::unnest(df[, c('id', 'item')], cols="item")
    # so <- vapply(df2$id, function(x) {
    #         paste(nn[as.numeric(strsplit(x, '/')[[1]])], collapse = " & ")
    #     }, character(1)
    # )


    # set_info <- data.frame(
    #     intersect = so, 
    #     Description = df2$item
    # )
    

    so <- vapply(split(d$Cluster, d$Description), 
        FUN = paste, 
        FUN.VALUE = character(1),
        collapse = " & ")

    set_info <- data.frame(
        intersect = so, 
        Description = names(so)
    )

    d2 <- merge(d, set_info, by="Description")
    n <- levels(d2$Cluster)
    cc <- yulab.utils::combinations(length(n))
    lv <- vapply(cc, function(i) paste(n[i], collapse = " & "), character(1))
    d2$intersect <- factor(d2$intersect, levels=lv)
    
    x@compareClusterResult <- d2

    return(x)
}

