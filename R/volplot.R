##' @rdname volplot
##' @exportMethod volplot
##' @author Guangchuang Yu
setMethod("volplot", signature(x = "enrichResult"),
        function(x, color = "zScore", 
                xintercept = 1, yintercept = 2,  
                showCategory = 5, label_format = 30,
                ...) {
            volplot.enrichResult(x = x, color = color, 
                xintercept = xintercept, yintercept = yintercept,
                showCategory = showCategory, label_format = label_format, 
                ...)
        })

##' @rdname volplot
##' @param font.size font size for `theme_dose()`
##' @param size font size to label selected categories specified by showCategory
volplot.enrichResult <- function(x, color = "zScore", 
                xintercept = 1, yintercept = 2, 
                showCategory=5, label_format = 30, 
                font.size=12, size = 5) {

    if (yintercept < 1) yintercept = -log10(yintercept)

    p <- ggplot(x@result, aes(x=log2(.data$FoldEnrichment), y= -log10(.data$p.adjust))) + 
        geom_point(aes(color=.data[[color]])) +
        geom_hline(yintercept = yintercept, lty='dashed') +
        geom_vline(xintercept = xintercept, lty='dashed')

    p <- p + set_enrichplot_color(type = "color", reverse = FALSE) +
        theme_dose(font.size) 

    if (is.numeric(showCategory)) {
        topN <- showCategory
        d <- dplyr::arrange(x@result, dplyr::desc(.data[[color]]))
        showCategory <- d$Description[1:topN]
    }

    label_func <- .label_format(label_format)
    p <- p + ggrepel::geom_text_repel(aes(label=label_func(.data$Description)), 
                    data = function(d) dplyr::filter(d, .data$Description %in% showCategory), 
                    size = size
        ) 
        
    p <- p + labs(x=bquote(paste(log[2], "(FoldEnrichment)")),
                y = bquote(paste(-log[10], "(p.adjust)"))
            )
    
    return(p)
}