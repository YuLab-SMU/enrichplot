##' upsetplot
##'
##'
##' @rdname upsetplot-methods
##' @aliases upsetplot,enrichResult,ANY-method
##' @param n number of categories to be plotted
##' @author Guangchuang Yu
##' @exportMethod upsetplot
##' @examples
##' require(DOSE)
##' data(geneList)
##' de=names(geneList)[1:100]
##' x <- enrichDO(de)
##' upsetplot(x, 8)
setMethod("upsetplot", signature(x="enrichResult"),
          function(x, n=10, ...) {
              upsetplot.enrichResult(x, n, ...)
          })

##' @rdname upsetplot-methods
##' @aliases upsetplot,gseaResult
##' @exportMethod upsetplot
setMethod("upsetplot", signature(x="gseaResult"),
          function(x, n=10, ...) {
              upsetplot.gseaResult(x, n, ...)
          })


##' @importFrom rlang check_installed
upsetplot.enrichResult <- function(x, n=10, ...) {
    df <- as.data.frame(x)
    id <- df$ID[1:n]
    des <- df$Description[1:n]
    glist <- geneInCategory(x)[id]
    names(glist) <- des
    ## g <- unique(unlist(glist))


    ## dat <- matrix(0, nrow=length(g), ncol=length(id))
    ## rownames(dat) <- g
    ## for (i in 1:length(id)) {
    ##     dat[glist[[i]], i] <- 1
    ## }
    ## colnames(dat) <- des

    ## ## cols <- ggtree:::color_scale("red", "blue")
    ## ## pv <- df$pvalue[1:n]
    ## ## idx <- sapply(pv, function(p) DOSE:::getIdx(p, min(pv), max(pv)))

    ## ## sets.bar.color = cols[idx],

    ## ## UpSetR <- "UpSetR"
    ## ## require(UpSetR, character.only = TRUE)
    ## ## upset <- eval(parse(text="upset"))

    ## upsetR::upset(as.data.frame(dat), nsets=n, ...)
    d <- list2df(glist)
    check_installed('tibble', 'for `upsetplot()`.')
    check_installed('ggupset', 'for `upsetplot()`.')
    res <- tibble::tibble(Description = split(d[,1], d[,2]))
    ggplot(res, aes_(x = ~Description)) + geom_bar() +
        theme_dose(font.size = 12) +
	xlab(NULL) + ylab(NULL) +
	ggupset::scale_x_upset(order_by = "freq")
}

##' @importFrom ggplot2 geom_violin
##' @importFrom ggplot2 geom_jitter
##' @importFrom rlang check_installed
upsetplot.gseaResult <- function(x, n = 10, type = "boxplot", ...) {
    n <- update_n(x, n)
    geneSets <- extract_geneSets(x, n)

    ## foldChange <- fc_readable(x, x@geneList)
    d <- list2df(geneSets)

    category <- split(d[,1], d[, 2])
    check_installed('tibble', 'for `upsetplot()`.')
    y <- tibble::tibble(Description = category,
                      gene = names(category),
                      foldChange = x@geneList[names(category)])

    if (type == "boxplot") {
        ly_dist <- geom_boxplot()
    } else {
        ly_dist <- geom_violin()
    }
    
    check_installed('ggupset', 'for `upsetplot()`.')
    ggplot(y, aes_(x = ~Description, y = ~foldChange)) +
        ly_dist +
        geom_jitter(width = .2, alpha = .6) +
        theme_dose(font.size = 12) +
        xlab(NULL) + ylab(NULL) +
        ggupset::scale_x_upset(order_by = "degree")
}

## @rdname upsetplot-methods
## @aliases upsetplot,compareClusterResult
## @exportMethod upsetplot
#setMethod("upsetplot", signature(x="compareClusterResult"),
#          function(x, n=10, ...) {
#              upsetplot.compareClusterResult(x, n, ...)
#          })


upsetplot.compareClusterResult <- function(x, n, ...) {
    x <- append_intersect(x)

    ## ggplot(x, aes(-10*log10(p.adjust), Description)) + geom_point() + facet_grid(set~., scales="free")

    ggplot(x, aes(.data$Cluster, .data$Description), showCategory=n) + 
        geom_point(aes(size=-10*log10(.data$p.adjust), color=.data$Cluster)) + 
        facet_grid(intersect ~ ., scales = "free", space = 'free') + guides(color = "none") +
        theme_dose(font.size = 12) +
        theme(strip.text = element_text(size = 14)) +
        xlab(NULL) + ylab(NULL) 
}

