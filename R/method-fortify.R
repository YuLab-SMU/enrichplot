##' convert enrichResult object for ggplot2
##'
##'
##' @title fortify
##' @param model enrichResult object
##' @param data not use here
##' @param showCategory Category numbers to show
##' @param by one of Count and GeneRatio
##' @param order logical
##' @param drop logical
##' @param split separate result by 'split' variable
##' @param ... additional parameter
##' @importFrom ggplot2 fortify
##' @method fortify enrichResult
##' @export
fortify.enrichResult <- function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, split=NULL, ...) {
    fortify.internal(model, data, showCategory, by, order, drop, split, ...)
}

##' @method fortify enrichResult
##' @export
fortify.gseaResult <- function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, split=NULL, ...) {
    fortify.internal(model, data, showCategory, by, order, drop, split, ...)
}


fortify.internal <- function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, split=NULL, ...) {
    res <- as.data.frame(model)
    if (inherits(model, "gseaResult")) {
        res$Count <- str_count(res$core_enrichment, "/") + 1
        res$.sign <- "activated"
        res$.sign[res$NES < 0] <- "suppressed"
    }
    if (drop) {
        res <- res[res$Count != 0, ]
    }
    if (inherits(model, "gseaResult")) {
        res$GeneRatio <- res$Count / res$setSize
    } else if (inherits(model, "enrichResult")) {
        res$GeneRatio <- parse_ratio(res$GeneRatio)
    }

    if (order) {
        if (by == "Count") {
            idx <- order(res$Count, decreasing=TRUE)
        } else {
            idx <- order(res$GeneRatio, decreasing=TRUE)
        }
        res <- res[idx,]
    }

    topN <- function(res, showCategory) {
        if ( is.numeric(showCategory) ) {
            if ( showCategory <= nrow(res) ) {
                res <- res[1:showCategory,]
            }
        } else { ## selected categories
            res <- res[res$ID %in% showCategory,]
        }
        return(res)
    }

    if (is.null(split)) {
        res <- topN(res, showCategory)
    } else {
        lres <- split(res, as.character(res[, split]))
        lres <- lapply(lres, topN, showCategory = showCategory)
        res <- do.call('rbind', lres)
    }

    res$Description <- factor(res$Description,
                              levels=rev(unique(res$Description)))

    return(res)
}

str_count <- function(string, pattern="") {
    sapply(string, str_count_item, pattern=pattern)
}

str_count_item <- function(string, pattern = "") {
    length(gregexpr(pattern, string)[[1]])
}

parse_ratio <- function(ratio) {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
    return(gsize/gcsize)
}


