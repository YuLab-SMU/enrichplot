##' upsetplot
##'
##'
##' @rdname upsetplot
##' @aliases upsetplot,enrichResult,ANY-method
##' @param n A number or a list of terms. If it is a number, 
##' the first n terms will be displayed. If it is a list of terms, 
##' the selected terms will be displayed.
##' @param enrichment_only whether only using enriched gene(for enrichResult) or core_enriched genes(for gseaResult) 
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

##' @rdname upsetplot
##' @aliases upsetplot,gseaResult
##' @exportMethod upsetplot
setMethod("upsetplot", signature(x="gseaResult"),
          function(x, n=10, ...) {
              upsetplot.gseaResult(x, n, ...)
          })

##' @rdname upsetplot
upsetplot.enrichResult <- function(x, n=10, enrichment_only = TRUE, ...) {
    # df <- as.data.frame(x)
    n <- update_n(x, n)
    glist <- get_geneSets(x = x, n = n, enrichment_only = enrichment_only)
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
    res <- tibble::tibble(Description = split(d[,1], d[,2]))
    ggplot(res, aes_(x = ~Description)) + geom_bar() +
        theme_dose(font.size = 12) +
        xlab(NULL) + ylab(NULL) +
        ggupset::scale_x_upset(order_by = "freq")
}

##' @rdname upsetplot
##' @importFrom ggplot2 geom_violin
##' @importFrom ggplot2 geom_jitter
##' @param type one of "boxplot" and "geom_violin"
##' @param nodeSize Size of nodes
upsetplot.gseaResult <- function(x, n = 10, type = "boxplot", nodeSize = 2, enrichment_only = TRUE, ...) {
    n <- update_n(x, n)
    geneSets <- get_geneSets(x = x, n = n, enrichment_only = enrichment_only)
    d <- list2df(geneSets)
    category <- split(d[,1], d[, 2])
    y <- tibble::tibble(Description = category,
          gene = names(category),
          foldChange = x@geneList[names(category)])


    if (type == "boxplot") {
        ly_dist <- geom_boxplot(outlier.size = nodeSize, ...)
    } else {
        ly_dist <- geom_violin(...)
    }

    ggplot(y, aes_(x = ~Description, y = ~foldChange)) +
        ly_dist +
        geom_jitter(width = .2, alpha = .6, size = nodeSize) +
        theme_dose(font.size = 12) +
        xlab(NULL) + ylab(NULL) +
        ggupset::scale_x_upset(order_by = "degree")
}


##' Get geneSets
##'
##' @param x enrichment result
##' @param n A number or a list of terms. If it is a number, 
##' the first n terms will be displayed. If it is a list of terms, 
##' the selected terms will be displayed.
##' @param enrichment_only whether only using enriched gene(for enrichResult) or core_enriched genes(for gseaResult) 
##' @noRd
get_geneSets <- function(x, n, enrichment_only) {
    if (enrichment_only) {
        geneSets <- extract_geneSets(x, n)
        ## foldChange <- fc_readable(x, x@geneList)
    } else {
        if (is.numeric(n)) {
            geneSets <- x@geneSets[x$ID[seq_len(n)]]
        } else {
            i <- match(n, x$Description)
            i <- i[!is.na(i)]
            nn <- x@result[i, "ID"]
            geneSets <- x@geneSets[nn]
            names(geneSets) <- x@result[i, "Description"]
        }
    }
    return(geneSets)
}









