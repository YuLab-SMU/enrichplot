##' @rdname ridgeplot
##' @exportMethod ridgeplot
setMethod("ridgeplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, fill = "p.adjust",
                   core_enrichment = TRUE, label_format = 30) {
              ridgeplot.gseaResult(x, showCategory = showCategory,
                                   fill = fill, core_enrichment = core_enrichment,
                                   label_format = label_format)
          })


##' @rdname ridgeplot
##' @importFrom ggplot2 scale_fill_gradientn
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 scale_x_reverse
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 scale_y_discrete
##' @author Guangchuang Yu
ridgeplot.gseaResult <- function(x, showCategory=30, fill="p.adjust",
                                 core_enrichment = TRUE, label_format = 30) {
    # has_package("ggridges")
    if (!is(x, "gseaResult"))
        stop("currently only support gseaResult")

    ## fill <- match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
    if (fill == "qvalue") {
        fill <- "qvalues"
    }
    if (!fill %in% colnames(x@result)) {
        stop("'fill' variable not available ...")
    }

    ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')

    n <- showCategory
    if (core_enrichment) {
        gs2id <- geneInCategory(x)[seq_len(n)]
    } else {
        gs2id <- x@geneSets[x$ID[seq_len(n)]]
    }

    gs2val <- lapply(gs2id, function(id) {
        res <- x@geneList[id]
        res <- res[!is.na(res)]
    })

    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]

    j <- order(x$NES[i], decreasing=FALSE)

    len <- sapply(gs2val, length)
    gs2val.df <- data.frame(category = rep(nn, times=len),
                            color = rep(x[i, fill], times=len),
                            value = unlist(gs2val))

    colnames(gs2val.df)[2] <- fill
    gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])

    label_func <- default_labeller(label_format)
    if(is.function(label_format)) {
        label_func <- label_format
    }

    ggplot(gs2val.df, aes_string(x="value", y="category", fill=fill)) +
        ggridges::geom_density_ridges() +
        ## scale_x_reverse() +
        scale_fill_continuous(low="red", high="blue", name = fill,
            guide=guide_colorbar(reverse=TRUE)) +
        scale_y_discrete(labels = label_func) +
        ## scale_fill_gradientn(name = fill, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        ## geom_vline(xintercept=0, color='firebrick', linetype='dashed') +
        xlab(NULL) + ylab(NULL) +  theme_dose()
}

