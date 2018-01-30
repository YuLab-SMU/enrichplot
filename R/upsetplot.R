
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

##' @importFrom UpSetR upset
upsetplot.enrichResult <- function(x, n=10, ...) {
    df <- as.data.frame(x)
    id <- df$ID[1:n]
    des <- df$Description[1:n]
    glist <- geneInCategory(x)[id]
    g <- unique(unlist(glist))


    dat <- matrix(0, nrow=length(g), ncol=length(id))
    rownames(dat) <- g
    for (i in 1:length(id)) {
        dat[glist[[i]], i] <- 1
    }
    colnames(dat) <- des

    ## cols <- ggtree:::color_scale("red", "blue")
    ## pv <- df$pvalue[1:n]
    ## idx <- sapply(pv, function(p) DOSE:::getIdx(p, min(pv), max(pv)))

    ## sets.bar.color = cols[idx],

    ## UpSetR <- "UpSetR"
    ## require(UpSetR, character.only = TRUE)
    ## upset <- eval(parse(text="upset"))

    upset(as.data.frame(dat), nsets=n, ...)
}

