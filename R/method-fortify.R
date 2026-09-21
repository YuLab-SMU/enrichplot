#' convert compareClusterResult to a data.frame that ready for plot
#'
#'
#' @rdname fortify
#' @title fortify
#' @param includeAll logical
#' @return data.frame
#' @importFrom ggplot2 fortify
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_head
#' @importFrom dplyr ungroup
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @export
#' @author Guangchuang Yu
fortify.compareClusterResult <- function(
    model,
    data,
    showCategory = 5,
    by = "geneRatio",
    split = NULL,
    includeAll = TRUE,
    ...
) {
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
    } else if (is.numeric(showCategory)) {
        topN <- function(res, showCategory) {
            if ("pvalue" %in% colnames(res)) {
                res <- arrange(res, .data$pvalue)
            } else {
                ## for groupGO
                res <- arrange(res, desc(.data$Count))
            }

            res %>% 
                group_by(.data$Cluster) %>% 
                slice_head(n = showCategory) %>% 
                ungroup() %>%
                as.data.frame()
        }

        if (!is.null(.split) && .split %in% colnames(clProf.df)) {
            lres <- split(clProf.df, as.character(clProf.df[, .split]))
            lres <- lapply(lres, topN, showCategory = showCategory)
            result <- as.data.frame(bind_rows(lres))
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
    GOlevel <- result[, c("ID", "Description")] ## GO ID and Term
    GOlevel <- unique(GOlevel)

    result <- result[result$Count != 0, ]
    result$Description <- factor(
        result$Description,
        levels = unique(rev(GOlevel[, 2]))
    )
    if (by == "rowPercentage") {
        Description <- Count <- NULL # to satisfy codetools
        result <- result %>%
            group_by(.data$Description) %>%
            mutate(
                Percentage = .data$Count / sum(.data$Count),
                Total = sum(.data$Count)
            ) %>%
            ungroup() %>%
            as.data.frame()

        ## label GO Description with gene counts.
        result$Description <- paste0(result$Description, " (", result$Total, ")")

        ## restore the original order of GO Description
        xx <- result[, c(2, 3)]
        xx <- unique(xx)
        rownames(xx) <- xx[, 1]
        Termlevel <- xx[as.character(GOlevel[, 1]), 2]

        ##drop the *Total* column
        result <- result[, colnames(result) != "Total"]

        result$Description <- factor(
            result$Description,
            levels = unique(rev(Termlevel))
        )
    } else if (by == "count") {
        result$GeneRatio <- yulab.utils::parse_ratio(result$GeneRatio)
    } else if (by == "geneRatio") {
        ## for result of ORA
        if (
            inherits(result$GeneRatio, "character") &&
                isTRUE(grepl("/", result$GeneRatio[1]))
        ) {
            gcsize <- as.numeric(sub(
                "^\\d+/",
                "",
                as.character(result$GeneRatio)
            ))
            result$GeneRatio <- yulab.utils::parse_ratio(result$GeneRatio)
            if (
                ("ONTOLOGY" %in% colnames(result)) &&
                    (length(unique(result$ONTOLOGY)) > 1)
            ) {
                # do nothing
            } else {
                cluster <- paste(
                    as.character(result$Cluster),
                    "\n",
                    "(",
                    gcsize,
                    ")",
                    sep = ""
                )
                orig_cls <- unique(result$Cluster)
                num_cls <- suppressWarnings(as.numeric(as.character(orig_cls)))
                
                if (any(is.na(num_cls))) {
                    idx <- order(orig_cls)
                } else {
                    idx <- order(num_cls)
                }
                
                lv <- unique(cluster)[idx]
                result$Cluster <- factor(cluster, levels = lv)
            }
        }
    } else {
        ## nothing
    }
    return(result)
}


#' convert enrichResult object for ggplot2
#'
#'
#' @title fortify
#' @rdname fortify
#' @param model 'enrichResult' or 'compareClusterResult' object
#' @param data not use here
#' @param showCategory Category numbers to show
#' @param by one of Count and GeneRatio
#' @param order logical
#' @param drop logical
#' @param split separate result by 'split' variable
#' @param ... additional parameter
#' @return data.frame
#' @importClassesFrom enrichit mnseaResult nseaResult
#' @importFrom ggplot2 fortify
## @method fortify enrichResult
#' @export
fortify.enrichResult <- function(
    model,
    data,
    showCategory = 5,
    by = "Count",
    order = FALSE,
    drop = FALSE,
    split = NULL,
    ...
) {
    fortify_internal(
        model = model,
        data = data,
        showCategory = showCategory,
        by = by,
        order = order,
        drop = drop,
        split = split
    )
}

## @method fortify gseaResult
#' @export
fortify.gseaResult <- function(
    model,
    data,
    showCategory = 5,
    by = "Count",
    order = FALSE,
    drop = FALSE,
    split = NULL,
    ...
) {
    fortify_internal(model, data, showCategory, by, order, drop, split)
}

#' @rdname fortify
## @method fortify nseaResult
#' @export
fortify.nseaResult <- function(
    model,
    data,
    showCategory = 5,
    by = "Count",
    order = FALSE,
    drop = FALSE,
    split = NULL,
    ...
) {
    fortify_internal(model, data, showCategory, by, order, drop, split)
}

#' @rdname fortify
#' @param level One of `"result"` or `"pathway"` for `mnseaResult`.
#' @param layer Optional layer or layers to retain for pathway-level output.
## @method fortify mnseaResult
#' @export
fortify.mnseaResult <- function(
    model,
    data,
    showCategory = 5,
    by = c("p.adjust", "NES", "contribution", "share"),
    order = FALSE,
    drop = FALSE,
    split = NULL,
    level = c("result", "pathway"),
    layer = NULL,
    ...
) {
    level <- match.arg(level)
    if (level == "result") {
        return(
            fortify_internal(
                model = model,
                data = data,
                showCategory = showCategory,
                by = "Count",
                order = order,
                drop = drop,
                split = split
            )
        )
    }

    by <- match.arg(by)
    df <- fortify_mnsea_contribution(model, level = "pathway")
    result_df <- .result_data(model)
    keep_cols <- intersect(
        c("ID", "Description", "pvalue", "p.adjust", "qvalue", "NES"),
        colnames(result_df)
    )

    if (length(keep_cols) > 0) {
        term_df <- unique(result_df[, keep_cols, drop = FALSE])
        join_by <- intersect(c("ID", "Description"), colnames(term_df))
        df <- merge(df, term_df, by = join_by, all.x = TRUE, sort = FALSE)
    }

    if (!is.null(layer)) {
        df <- df[df$layer %in% layer, , drop = FALSE]
    }

    if (nrow(df) == 0) {
        return(df)
    }

    if ("NES" %in% colnames(df)) {
        df$.sign <- "activated"
        df$.sign[df$NES < 0] <- "suppressed"
    }

    if (is.numeric(showCategory)) {
        rank_df <- .rank_mnsea_terms(df, by = by)
        show_n <- min(showCategory, nrow(rank_df))
        rank_df <- rank_df[seq_len(show_n), , drop = FALSE]
        df <- df[df$ID %in% rank_df$ID, , drop = FALSE]
        term_levels <- rank_df$Description
    } else {
        keep <- unique(showCategory)
        df <- df[df$Description %in% keep | df$ID %in% keep, , drop = FALSE]
        term_levels <- unique(df$Description)
    }

    if (by %in% c("contribution", "share")) {
        idx <- order(df[[by]], decreasing = TRUE)
        df <- df[idx, , drop = FALSE]
    }

    df$Description <- factor(df$Description, levels = rev(unique(term_levels)))
    rownames(df) <- NULL
    df
}


fortify_internal <- function(
    model,
    data,
    showCategory = 5,
    by = "Count",
    order = FALSE,
    drop = FALSE,
    split = NULL
) {
    res <- .result_data(model)
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
            idx <- order(res$Count, decreasing = TRUE)
        } else {
            idx <- order(res$GeneRatio, decreasing = TRUE)
        }
        res <- res[idx, ]
    }

    topN <- function(res, showCategory) {
        if (is.numeric(showCategory)) {
            if (showCategory <= nrow(res)) {
                res <- res[1:showCategory, ]
            }
        } else {
            ## selected categories
            res <- res[res$Description %in% showCategory, ]
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

    res$Description <- factor(
        res$Description,
        levels = rev(unique(res$Description))
    )

    return(res)
}

str_count <- function(string, pattern = "") {
    sapply(string, FUN = function(i) {
        length(unlist(strsplit(i, split = pattern)))
    })
}

parse_ratio <- function(ratio) {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
    return(gsize / gcsize)
}
