#' Upsetplot
#'
#' Upsetplot
#'
#' @rdname upsetplot-methods
#' @aliases upsetplot,enrichResult,ANY-method
#' @param n number of categories to be plotted
#' @param ... additional parameters
#' @author Guangchuang Yu
#' @exportMethod upsetplot
#' @examples
#' \dontrun{
#' library(DOSE)
#' data(geneList)
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' upsetplot(x, 8)
#' }
setMethod("upsetplot", signature(x="enrichResult"),
          function(
              x,
              n = 10,
              type = "boxplot",
              layer = NULL,
              value = c("score", "abs_score"),
              core_enrichment = FALSE,
              ...
          ) {
              upsetplot.enrichResult(x, n, ...)
          })

#' @rdname upsetplot-methods
#' @aliases upsetplot,gseaResult
#' @exportMethod upsetplot
setMethod("upsetplot", signature(x="gseaResult"),
          function(
              x,
              n = 10,
              type = "boxplot",
              layer = NULL,
              value = c("score", "abs_score"),
              core_enrichment = FALSE,
              ...
          ) {
              upsetplot.gseaResult(x, n, type = type, ...)
          })

#' @rdname upsetplot-methods
#' @aliases upsetplot,mnseaResult
#' @param layer Optional `mnsea` layer. When `NULL`, use collapsed scores.
#' @param value score summary to display for overlapping features.
#' @param core_enrichment logical. Should only core mnsea features be used?
#' @exportMethod upsetplot
setMethod("upsetplot", signature(x="mnseaResult"),
          function(
              x,
              n = 10,
              type = "boxplot",
              layer = NULL,
              value = c("score", "abs_score"),
              core_enrichment = FALSE,
              ...
          ) {
              upsetplot.mnseaResult(
                  x,
                  n = n,
                  type = type,
                  layer = layer,
                  value = value,
                  core_enrichment = core_enrichment,
                  ...
              )
          })


#' @importFrom rlang check_installed
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
    require_suggested(c('tibble', 'ggupset'), 'for `upsetplot()`.')
    res <- tibble::tibble(Description = split(d[,1], d[,2]))
    ggplot(res, aes(x = .data$Description)) + geom_bar() +
        theme_dose(font.size = 12) +
	xlab(NULL) + ylab(NULL) +
	ggupset::scale_x_upset(order_by = "freq")
}

#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 geom_jitter
#' @importFrom rlang check_installed
upsetplot.gseaResult <- function(x, n = 10, type = "boxplot", ...) {
    n <- update_n(x, n)
    geneSets <- extract_geneSets(x, n)
    labels <- get_geneSet_labels(geneSets)

    ## foldChange <- fc_readable(x, x@geneList)
    d <- list2df(geneSets)
    d$Description <- labels[as.character(d$categoryID)]

    category <- split(d$Description, d$Gene)
    require_suggested('tibble', 'for `upsetplot()`.')
    y <- tibble::tibble(Description = category,
                      gene = names(category),
                      foldChange = x@geneList[names(category)])

    if (type == "boxplot") {
        ly_dist <- geom_boxplot()
    } else {
        ly_dist <- geom_violin()
    }
    
    require_suggested('ggupset', 'for `upsetplot()`.')
    ggplot(y, aes(x = .data$Description, y = .data$foldChange)) +
        ly_dist +
        geom_jitter(width = .2, alpha = .6) +
        theme_dose(font.size = 12) +
        xlab(NULL) + ylab(NULL) +
        ggupset::scale_x_upset(order_by = "degree")
}

upsetplot.mnseaResult <- function(
    x,
    n = 10,
    type = "boxplot",
    layer = NULL,
    value = c("score", "abs_score"),
    core_enrichment = FALSE,
    ...
) {
    value <- match.arg(value)
    y <- build_mnsea_upset_df(
        x,
        n = n,
        layer = layer,
        value = value,
        core_enrichment = core_enrichment
    )

    if (type == "boxplot") {
        ly_dist <- geom_boxplot()
    } else {
        ly_dist <- geom_violin()
    }

    require_suggested("ggupset", "for `upsetplot()`.")
    ggplot(y, aes(x = .data$Description, y = .data[[value]])) +
        ly_dist +
        geom_jitter(width = .2, alpha = .6) +
        theme_dose(font.size = 12) +
        xlab(NULL) + ylab(mnsea_plot_label(value)) +
        ggupset::scale_x_upset(order_by = "degree")
}

build_mnsea_upset_df <- function(
    object,
    n = 10,
    layer = NULL,
    value = c("score", "abs_score"),
    core_enrichment = FALSE
) {
    value <- match.arg(value)
    selected_ids <- resolve_mnsea_ridge_ids(object, showCategory = n)
    geneList <- get_mnsea_ranked_scores(object, layer = layer)

    feature_df <- do.call(
        rbind,
        lapply(selected_ids, function(id) {
            df <- fortify_mnsea_contribution(
                object,
                level = "feature",
                pathway_id = id,
                layer = layer
            )
            if (core_enrichment && "is_core" %in% colnames(df)) {
                df <- df[df$is_core %in% TRUE, , drop = FALSE]
            }
            df
        })
    )

    if (is.null(feature_df) || nrow(feature_df) == 0) {
        yulab.utils::yulab_abort(
            "No mnsea features available for upsetplot after filtering."
        )
    }

    feature_df <- feature_df[feature_df$ID %in% selected_ids, , drop = FALSE]
    feature_df$Feature <- as.character(feature_df$Feature)
    feature_df <- feature_df[feature_df$Feature %in% names(geneList), , drop = FALSE]
    if (nrow(feature_df) == 0) {
        yulab.utils::yulab_abort(
            "No overlap between mnsea features and ranked scores for upsetplot."
        )
    }

    labels <- get_term_labels(object, selected_ids)
    feature_df$Description <- unname(as.character(labels[feature_df$ID]))
    feature_df <- feature_df[!duplicated(feature_df[, c("Feature", "ID")]), , drop = FALSE]

    category <- split(feature_df$Description, feature_df$Feature)
    feature_ids <- names(category)
    score <- geneList[feature_ids]
    y <- data.frame(
        Feature = feature_ids,
        score = as.numeric(score),
        abs_score = abs(as.numeric(score)),
        layer = if (is.null(layer)) "collapsed" else as.character(layer),
        stringsAsFactors = FALSE
    )
    y$Description <- unname(category[y$Feature])
    y
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

    ggplot(x, aes(x = .data$Cluster, y = .data$Description), showCategory=n) + 
        geom_point(aes(size = -10 * log10(.data$p.adjust), color = .data$Cluster)) + 
        facet_grid(intersect ~ ., scales = "free", space = 'free') + guides(color = "none") +
        theme_dose(font.size = 12) +
        theme(strip.text = element_text(size = 14)) +
        xlab(NULL) + ylab(NULL) 
}

