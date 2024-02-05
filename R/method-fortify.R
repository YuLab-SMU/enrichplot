

##' convert compareClusterResult to a data.frame that ready for plot
##'
##'
##' @rdname fortify
##' @title fortify
##' @param includeAll logical
##' @return data.frame
##' @importFrom ggplot2 fortify
##' @importFrom plyr ddply
##' @importFrom plyr mdply
##' @importFrom plyr .
##' @method fortify compareClusterResult
##' @export
##' @author Guangchuang Yu
fortify.compareClusterResult <- function(model, data, showCategory=5,
                                         by="geneRatio", split=NULL,
                                         includeAll=TRUE, ...) {
    clProf.df <- as.data.frame(model)
    .split <- split
    if ("core_enrichment" %in% colnames(clProf.df)) {
        clProf.df$Count <- str_count(clProf.df$core_enrichment, "/")
        clProf.df$.sign <- "activated"
        clProf.df$.sign[clProf.df$NES < 0] <- "suppressed"
        clProf.df$GeneRatio <- clProf.df$Count / clProf.df$setSize
    }
    ## get top 5 (default) categories of each gene cluster.
    if (is.null(showCategory)) {
        result <- clProf.df
    } else if(is.numeric(showCategory)){
        Cluster <- NULL # to satisfy codetools

        topN <- function(res, showCategory) {
            ddply(.data = res,
                  .variables = .(Cluster),
                  .fun = function(df, N) {
                      if (length(df$Count) > N) {
                          if (any(colnames(df) == "pvalue")) {
                              idx <- order(df$pvalue, decreasing=FALSE)[1:N]
                          } else {
                              ## for groupGO
                              idx <- order(df$Count, decreasing=T)[1:N]
                          }
                          return(df[idx,])
                      } else {
                          return(df)
                      }
                  },
                  N=showCategory
                  )

        }

        if (!is.null(.split) && .split %in% colnames(clProf.df)) {
            lres <- split(clProf.df, as.character(clProf.df[, .split]))
            lres <- lapply(lres, topN, showCategory = showCategory)
            result <- do.call('rbind', lres)
        } else {
            result <- topN(clProf.df, showCategory)
        }

    } else {
        result <- subset(clProf.df, Description %in% showCategory)
    }

    ID <- NULL
    if (includeAll == TRUE) {
        result <- subset(clProf.df, ID %in% result$ID)
    }

    ## remove zero count
    result$Description <- as.character(result$Description) ## un-factor
    GOlevel <- result[,c("ID", "Description")] ## GO ID and Term
    GOlevel <- unique(GOlevel)



    result <- result[result$Count != 0, ]
    result$Description <- factor(result$Description,
                                 levels=rev(GOlevel[,2]))
    if (by=="rowPercentage") {
        Description <- Count <- NULL # to satisfy codetools
        result <- ddply(result,
                        .(Description),
                        transform,
                        Percentage = Count/sum(Count),
                        Total = sum(Count))

        ## label GO Description with gene counts.
        x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
        y <- sapply(x[,3], paste, ")", sep="")
        result$Description <- y

        ## restore the original order of GO Description
        xx <- result[,c(2,3)]
        xx <- unique(xx)
        rownames(xx) <- xx[,1]
        Termlevel <- xx[as.character(GOlevel[,1]),2]

        ##drop the *Total* column
        result <- result[, colnames(result) != "Total"]

        result$Description <- factor(result$Description,
                                     levels=rev(Termlevel))

    } else if (by == "count") {
        ## nothing
    } else if (by == "geneRatio") {
        ## for result of ORA
        # if (class(result$GeneRatio) == "character" && grep("/", result$GeneRatio[1])) {
        if (inherits(result$GeneRatio, "character") && grep("/", result$GeneRatio[1])) {
            gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
            gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
            result$GeneRatio <- gsize/gcsize
            if (("ONTOLOGY" %in% colnames(result)) && (length(unique(result$ONTOLOGY)) > 1)){
                # do nothing
            } else {
                cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")",
                                 sep="")
                lv <- unique(cluster)[order(as.numeric(unique(result$Cluster)))]
                result$Cluster <- factor(cluster, levels = lv)
            }
        }

    } else {
        ## nothing
    }
    return(result)
}


##' convert enrichResult object for ggplot2
##'
##'
##' @title fortify
##' @rdname fortify
##' @param model 'enrichResult' or 'compareClusterResult' object
##' @param data not use here
##' @param showCategory Category numbers to show
##' @param by one of Count and GeneRatio
##' @param order logical
##' @param drop logical
##' @param split separate result by 'split' variable
##' @param ... additional parameter
##' @return data.frame
##' @importFrom ggplot2 fortify
##' @method fortify enrichResult
##' @export
fortify.enrichResult <- function(model, data, showCategory=5, by = "Count",
                                 order=FALSE, drop=FALSE, split=NULL, ...) {
    fortify_internal(model, data, showCategory, by, order, drop, split, ...)
}

##' @method fortify gseaResult
##' @export
fortify.gseaResult <- function(model, data, showCategory=5, by = "Count",
                               order=FALSE, drop=FALSE, split=NULL, ...) {
    fortify_internal(model, data, showCategory, by, order, drop, split, ...)
}


fortify_internal <- function(model, data, showCategory=5, by = "Count",
                             order=FALSE, drop=FALSE, split=NULL, ...) {
    res <- as.data.frame(model)
    res <- res[!is.na(res$Description), ]
    if (inherits(model, "gseaResult")) {
        res$Count <- str_count(res$core_enrichment, "/")
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
        if ("BgRatio" %in% colnames(res)) {
            ## groupGO output doesn't have this column
            res$BgRatio <- parse_ratio(res$BgRatio)
        }
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
            res <- res[res$Description %in% showCategory,]
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
    sapply(string, FUN = function(i) length(unlist(strsplit(i, split = pattern))))
}

parse_ratio <- function(ratio) {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
    return(gsize/gcsize)
}


