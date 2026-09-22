#' @rdname ridgeplot
#' @exportMethod ridgeplot
setMethod(
    "ridgeplot",
    signature(x = "gseaResult"),
    function(
        x,
        showCategory = 30,
        fill = "p.adjust",
        core_enrichment = TRUE,
        label_format = 30,
        ...
    ) {
        ridgeplot.gseaResult(
            x,
            showCategory = showCategory,
            fill = fill,
            core_enrichment = core_enrichment,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname ridgeplot
#' @param layer Optional `mnsea` layer. When `NULL`, use collapsed scores.
#' @exportMethod ridgeplot
setMethod(
    "ridgeplot",
    signature(x = "mnseaResult"),
    function(
        x,
        showCategory = 30,
        fill = "p.adjust",
        core_enrichment = TRUE,
        label_format = 30,
        orderBy = "NES",
        decreasing = FALSE,
        stat = "density_ridges",
        layer = NULL,
        ...
    ) {
        ridgeplot.mnseaResult(
            x,
            showCategory = showCategory,
            fill = fill,
            core_enrichment = core_enrichment,
            label_format = label_format,
            orderBy = orderBy,
            decreasing = decreasing,
            stat = stat,
            layer = layer,
            ...
        )
    }
)


#' @rdname ridgeplot
#' @param orderBy The order of the Y-axis
#' @param decreasing logical. Should the orderBy order be increasing or decreasing?
#' @param stat statistic passed to `ggridges::geom_density_ridges()`.
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_x_reverse
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom rlang check_installed
#' @importFrom yulab.utils yulab_abort
#' @importFrom yulab.utils yulab_warn
#' @author Guangchuang Yu
ridgeplot.gseaResult <- function(
    x,
    showCategory = 30,
    fill = "p.adjust",
    core_enrichment = TRUE,
    label_format = 30,
    orderBy = "NES",
    decreasing = FALSE,
    stat = "density_ridges"
) {
    ## Input validation with better error messages
    check_input(x, type = "gseaResult", arg_name = "x")
    
    if (!fill %in% colnames(x@result)) {
        yulab_abort(paste0("'", fill, "' variable not available in result"), 
                        class = "missing_column_error")
    }

    ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')
    if (orderBy != 'NES' && !orderBy %in% colnames(x@result)) {
        yulab_warn('wrong orderBy parameter; set to default `orderBy = "NES"`',
                     class = "parameter_warning")
        orderBy <- "NES"
    }
    
    ## Optimized category selection
    if (inherits(showCategory, 'numeric')) {
        selected <- seq_len(min(showCategory, nrow(x@result)))
    } else if (inherits(showCategory, "character")) {
        ii <- match(showCategory, x@result$Description)
        if (all(is.na(ii))) {
            ii <- match(showCategory, x@result$ID)
        }
        ii <- ii[!is.na(ii)]
        if (length(ii) == 0) {
            yulab_warn("No matching categories found, using first 10",
                          class = "category_warning")
            ii <- seq_len(min(10, nrow(x@result)))
        }
        selected <- x@result[ii, "ID"]
    } else {
        yulab_warn("showCategory should be a number of pathways or a vector of selected pathways",
                       class = "parameter_warning")
        selected <- seq_len(min(10, nrow(x@result)))
    }

    ## Optimized gene set extraction
    if (core_enrichment) {
        gs2id <- geneInCategory(x)[selected]
    } else {
        gs2id <- x@geneSets[names(x@geneSets) %in% selected]
    }

    ## Optimized gene name mapping
    if (x@readable && length(x@gene2Symbol) > 0) {
        gene_names <- names(x@geneList)
        symbol_match <- match(gene_names, names(x@gene2Symbol))
        valid_matches <- !is.na(symbol_match)
        names(x@geneList)[valid_matches] <- x@gene2Symbol[symbol_match[valid_matches]]
    }

    ## Vectorized data preparation
    gs2val <- lapply(gs2id, function(id) {
        res <- x@geneList[id]
        res[!is.na(res)]
    })

    ids <- names(gs2val)
    i <- match(ids, x$ID)
    labels <- x$Description[i]
    fill_values <- x[i, fill]
    order_values <- x@result[[orderBy]][i]

    pruned <- prune_ridgeplot_groups(
        gs2val = gs2val,
        ids = ids,
        labels = labels,
        fill_values = fill_values,
        order_values = order_values
    )
    gs2val <- pruned$gs2val
    ids <- pruned$ids
    labels <- pruned$labels
    fill_values <- pruned$fill_values
    order_values <- pruned$order_values

    ## Optimized ordering
    j <- order(order_values, decreasing = decreasing)
    
    ## Efficient data frame construction
    len <- lengths(gs2val)
    
    gs2val.df <- data.frame(
        category = rep(labels, times = len),
        color = rep(fill_values, times = len),
        value = unlist(gs2val, use.names = FALSE)
    )

    colnames(gs2val.df)[2] <- fill
    gs2val.df$category <- factor(gs2val.df$category, levels = labels[j])

    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    require_suggested('ggridges', 'for `ridgeplot()`.')

    ggplot(
        gs2val.df,
        aes(x = .data[["value"]], y = .data[["category"]], fill = .data[[fill]])
    ) +
        ggridges::geom_density_ridges(stat = stat) +
        set_enrichplot_color(type = "fill", name = fill, transform = 'log10') +
        scale_y_discrete(labels = label_func) +
        xlab(NULL) +
        ylab(NULL) +
        theme_dose()
}

ridgeplot.mnseaResult <- function(
    x,
    showCategory = 30,
    fill = "p.adjust",
    core_enrichment = TRUE,
    label_format = 30,
    orderBy = "NES",
    decreasing = FALSE,
    stat = "density_ridges",
    layer = NULL
) {
    if (!fill %in% colnames(x@result)) {
        yulab_abort(
            paste0("'", fill, "' variable not available in result"),
            class = "missing_column_error"
        )
    }

    if (orderBy != "NES" && !orderBy %in% colnames(x@result)) {
        yulab_warn(
            'wrong orderBy parameter; set to default `orderBy = "NES"`',
            class = "parameter_warning"
        )
        orderBy <- "NES"
    }

    gs2val.df <- build_mnsea_ridge_df(
        x,
        showCategory = showCategory,
        fill = fill,
        core_enrichment = core_enrichment,
        orderBy = orderBy,
        decreasing = decreasing,
        layer = layer
    )

    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    require_suggested("ggridges", "for `ridgeplot()`.")

    ggplot(
        gs2val.df,
        aes(x = .data[["value"]], y = .data[["category"]], fill = .data[[fill]])
    ) +
        ggridges::geom_density_ridges(stat = stat) +
        set_enrichplot_color(type = "fill", name = fill, transform = "log10") +
        scale_y_discrete(labels = label_func) +
        xlab(NULL) +
        ylab(NULL) +
        theme_dose()
}

build_mnsea_ridge_df <- function(
    object,
    showCategory = 30,
    fill = "p.adjust",
    core_enrichment = TRUE,
    orderBy = "NES",
    decreasing = FALSE,
    layer = NULL
) {
    result_df <- .result_data(object)
    if (nrow(result_df) == 0) {
        yulab_abort("No mnsea pathways available for `ridgeplot()`.")
    }

    selected_ids <- resolve_mnsea_ridge_ids(object, showCategory = showCategory)
    geneList <- get_mnsea_ranked_scores(object, layer = layer)

    feature_list <- lapply(selected_ids, function(id) {
        feature_df <- fortify_mnsea_contribution(
            object,
            level = "feature",
            pathway_id = id,
            layer = layer
        )
        if (core_enrichment && "is_core" %in% colnames(feature_df)) {
            feature_df <- feature_df[isTRUE(feature_df$is_core) | feature_df$is_core %in% TRUE, , drop = FALSE]
        }
        feature_df
    })
    names(feature_list) <- selected_ids

    gs2val <- lapply(feature_list, function(feature_df) {
        feature_ids <- unique(as.character(feature_df$Feature))
        feature_ids <- feature_ids[!is.na(feature_ids) & nzchar(feature_ids)]
        feature_ids <- intersect(feature_ids, names(geneList))
        geneList[feature_ids]
    })

    keep <- lengths(gs2val) > 0
    if (!any(keep)) {
        yulab_abort(
            "No overlap between mnsea features and ranked scores for the selected pathway/layer."
        )
    }
    gs2val <- gs2val[keep]
    selected_ids <- names(gs2val)

    labels <- get_term_labels(object, selected_ids)
    label_vec <- unname(as.character(labels[selected_ids]))
    fill_values <- object@result[selected_ids, fill]
    order_values <- object@result[selected_ids, orderBy]

    ## mnsea mechanism plots are built from per-feature layer contributions,
    ## which are much smaller than full ranked lists. Only drop empty groups
    ## (already done above); do not impose the >= 3 points that gseaResult
    ## requires for density estimation.
    pruned <- prune_ridgeplot_groups(
        gs2val = gs2val,
        ids = selected_ids,
        labels = label_vec,
        fill_values = fill_values,
        order_values = order_values,
        min_size = 1L
    )
    gs2val <- pruned$gs2val
    selected_ids <- pruned$ids
    label_vec <- pruned$labels
    fill_values <- pruned$fill_values
    order_values <- pruned$order_values

    j <- order(order_values, decreasing = decreasing)

    len <- lengths(gs2val)
    gs2val.df <- data.frame(
        ID = rep(selected_ids, times = len),
        category = rep(label_vec, times = len),
        value = unlist(gs2val, use.names = FALSE),
        layer = if (is.null(layer)) "collapsed" else as.character(layer),
        stringsAsFactors = FALSE
    )
    gs2val.df[[fill]] <- rep(fill_values, times = len)
    gs2val.df$category <- factor(gs2val.df$category, levels = label_vec[j])
    gs2val.df
}

prune_ridgeplot_groups <- function(
    gs2val,
    ids,
    labels,
    fill_values,
    order_values,
    min_size = 3L
) {
    counts <- lengths(gs2val)
    keep <- counts >= min_size

    if (!all(keep)) {
        yulab_warn(
            paste0(
                "Dropping ",
                sum(!keep),
                " gene set(s) with fewer than ",
                min_size,
                " ranked values for `ridgeplot()`."
            ),
            class = "insufficient_values_warning"
        )
    }

    if (!any(keep)) {
        yulab_abort(
            paste0(
                "No selected gene sets have at least ",
                min_size,
                " ranked values for `ridgeplot()`."
            )
        )
    }

    list(
        gs2val = gs2val[keep],
        ids = ids[keep],
        labels = labels[keep],
        fill_values = fill_values[keep],
        order_values = order_values[keep]
    )
}

resolve_mnsea_ridge_ids <- function(object, showCategory = 30) {
    result_df <- .result_data(object)
    if (is.numeric(showCategory)) {
        n <- min(as.integer(showCategory), nrow(result_df))
        if (is.na(n) || n < 1) {
            yulab_warn(
                "showCategory should be a positive number of pathways or a vector of selected pathways",
                class = "parameter_warning"
            )
            n <- min(10, nrow(result_df))
        }
        return(as.character(result_df$ID[seq_len(n)]))
    }

    if (is.character(showCategory)) {
        idx <- resolve_term_rows(object, showCategory)
        if (length(idx) == 0) {
            yulab_warn(
                "No matching categories found, using first 10",
                class = "category_warning"
            )
            return(as.character(result_df$ID[seq_len(min(10, nrow(result_df)))]))
        }
        return(as.character(result_df$ID[idx]))
    }

    yulab_warn(
        "showCategory should be a number of pathways or a vector of selected pathways",
        class = "parameter_warning"
    )
    as.character(result_df$ID[seq_len(min(10, nrow(result_df)))])
}
