##' @rdname heatplot
##' @exportMethod heatplot
setMethod("heatplot", signature(x = "enrichResult"),
          function(x, showCategory = 30, ...) {
              heatplot.enrichResult(x, showCategory, ...)
          })

##' @rdname heatplot
##' @exportMethod heatplot
setMethod("heatplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, ...) {
              heatplot.enrichResult(x, showCategory, ...)
          })



##' @rdname heatplot
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 theme_minimal
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 scale_y_discrete
##' @importFrom ggplot2 scale_fill_gradient2
##' @importFrom rlang check_installed
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param symbol symbol of the nodes, one of "rect"(the default) and "dot"
##' by default wraps names longer that 30 characters
##' @param pvalue pvalue of genes
##' @author Guangchuang Yu
heatplot.enrichResult <- function(x, showCategory = 30, symbol = "rect", foldChange = NULL,
                                  pvalue = NULL, label_format = 30) {

    symbol <- match.arg(symbol, c("rect", "dot"))
    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }

    n <- update_n(x, showCategory)
    geneSets <- extract_geneSets(x, n)
    foldChange <- fc_readable(x, foldChange)
    pvalue <- fc_readable(x, pvalue)
    d <- list2df(geneSets)
    if (!is.null(foldChange)) {
        d$foldChange <- foldChange[as.character(d[,2])]
    } 

    if (!is.null(pvalue)) {
        d$pvalue <- pvalue[as.character(d[,2])]
    } 
 
    p <- ggplot(d, aes_(~Gene, ~categoryID))

    if (symbol == "rect") {
        p <- p + geom_tile(color = 'white')
    } 

    get_dotp <-function(p, foldChange, pvalue) {
        if (is.null(foldChange) & is.null(pvalue)) {
            p <- p + geom_point(color = 'black', shape = 21, fill = "black", size = 5)
            return(p)
        }

        if (!is.null(foldChange) & !is.null(pvalue)) {
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
   
        check_installed('scales', 'for `heatplot()`.')
    
        inv <- function(x) base^(-x)
        scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
                  scales::log_breaks(base = base), 
                  domain = c(1e-100, Inf))
    }

    if (symbol == "dot") {
        p <- get_dotp(p, foldChange, pvalue)     
        ## only dot need size(pvalue) parameter
        if (!is.null(pvalue)) {
            p <- p + aes_(size = ~pvalue) + 
                scale_size_continuous(range=c(3, 8), 
                    trans = reverselog_trans(10))
        } 
    }

    if (!is.null(foldChange)) {
            p <- p + aes_(fill = ~foldChange) + 
                set_enrichplot_color(colors = get_enrichplot_color(3), type = "fill")
                # scale_fill_gradient2(name = "fold change", low = "#327eba",
                #                    mid = "white", high = "#e06663") +

    }
    

    p + xlab(NULL) + ylab(NULL) + theme_minimal() +
        scale_y_discrete(labels = label_func) +
        theme(panel.grid.major = element_blank(),
              axis.text.x=element_text(angle = 60, hjust = 1))
}

