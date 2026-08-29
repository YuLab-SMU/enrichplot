#' @rdname pairwise_termsim
#' @exportMethod pairwise_termsim
setMethod("pairwise_termsim", signature(x = "enrichResult"),
    function(x, method = "JC", semData = NULL, showCategory = NULL, layer = NULL) {
        pairwise_termsim.enrichResult(x, method = method,
            semData = semData, showCategory = showCategory)
    })

#' @rdname pairwise_termsim
#' @exportMethod pairwise_termsim
setMethod("pairwise_termsim", signature(x = "gseaResult"),
    function(x, method = "JC", semData = NULL, showCategory = NULL, layer = NULL) {
        pairwise_termsim.enrichResult(x, method = method,
            semData = semData, showCategory = showCategory)
    })

#' @rdname pairwise_termsim
#' @exportMethod pairwise_termsim
setMethod("pairwise_termsim", signature(x = "mnseaResult"),
    function(x, method = "JC", semData = NULL, showCategory = NULL, layer = NULL) {
        if (is.null(showCategory)) {
            showCategory <- .default_pairwise_termsim_category(x)
        }
        n <- showCategory
        if (is.numeric(n) && n == 0) {
            stop("no enriched term found...")
        }
        if (!is.numeric(n) && length(n) == 0) {
            stop("no enriched term found...")
        }

        selected <- select_terms(x, showCategory)
        if (nrow(selected$result) == 0) {
            stop("no enriched term found...")
        }

        feature_df <- prepare_emapplot_mnsea_feature_data(
            x,
            ids = selected$ids,
            layer = layer
        )
        if (nrow(feature_df) == 0) {
            stop("no mnsea features available for pairwise termsim.")
        }

        geneSets <- lapply(selected$ids, function(id) {
            unique(as.character(feature_df$Feature[feature_df$ID == id]))
        })
        names(geneSets) <- selected$ids
        x@termsim <- get_similarity_matrix(
            y = selected$result,
            geneSets = geneSets,
            method = method,
            semData = semData
        )
        x@method <- method
        return(x)
    })

#' @rdname pairwise_termsim
#' @exportMethod pairwise_termsim
setMethod("pairwise_termsim", signature(x = "compareClusterResult"),
    function(x, method = "JC", semData = NULL, showCategory = NULL, layer = NULL) {
        pairwise_termsim.compareClusterResult(x, method = method,
            semData = semData, showCategory = showCategory)
    })


#' @rdname pairwise_termsim
prepare_pairwise_termsim_data <- function(x, showCategory) {
    selected <- select_terms(x, showCategory)
    list(
        result = selected$result,
        geneSets = geneInCategory(x),
        selection = selected$selection
    )
}

#' @rdname pairwise_termsim
pairwise_termsim.enrichResult <- function(x, method = "JC", semData = NULL, showCategory = NULL) {
    if (is.null(showCategory)) {
        showCategory <- .default_pairwise_termsim_category(x)
    }

    selected <- prepare_pairwise_termsim_data(x, showCategory)
    n <- selected$selection
    if (is.numeric(n) && n == 0) {
        stop("no enriched term found...")
    }
    if (!is.numeric(n) && length(n) == 0) {
        stop("no enriched term found...")
    }
    y <- selected$result
    geneSets <- selected$geneSets

    x@termsim <- get_similarity_matrix(y = y, geneSets = geneSets, method = method,
                semData = semData)
    x@method <- method
    return(x)
}


#' @rdname pairwise_termsim
pairwise_termsim.compareClusterResult <- function(x, method = "JC", semData = NULL, 
                                                  showCategory = NULL) {
    if (is.null(showCategory)) {
        showCategory <- .default_pairwise_termsim_category(x)
    }

    y <- fortify(x, showCategory=showCategory, includeAll=TRUE, split=NULL)
    y$Cluster <- sub("\n.*", "", y$Cluster)
    ## y_union <- get_y_union(y = y, showCategory = showCategory)
    if ("core_enrichment" %in% colnames(y)) {
        y$geneID <- y$core_enrichment
    }
    y_union <- merge_compareClusterResult(y)
    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/",
                                  fixed = TRUE), 
                         y_union$ID)
    x@termsim <- get_similarity_matrix(y = y_union, geneSets = geneSets, method = method,
                semData = semData)                              
    x@method <- method
    return(x)    
}


.default_pairwise_termsim_category <- function(x, min_default = 200) {
    min(nrow(x), min_default)
}
