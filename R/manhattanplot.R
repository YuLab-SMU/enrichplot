#' @rdname manhattanplot
#' @exportMethod manhattanplot
#' @author Guangchuang Yu
setMethod(
    "manhattanplot",
    signature(x = "enrichResult"),
    function(
        x,
        color = "p.adjust",
        showCategory = 5,
        size = "Count",
        split = NULL,
        font.size = 12,
        title = "",
        label_format = 30,
        ...
    ) {
        manhattanplot.enrichResult(
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname manhattanplot
#' @exportMethod manhattanplot
setMethod(
    "manhattanplot",
    signature(x = "gseaResult"),
    function(
        x,
        color = "p.adjust",
        showCategory = 5,
        size = "Count",
        split = NULL,
        font.size = 12,
        title = "",
        label_format = 30,
        ...
    ) {
        manhattanplot.enrichResult(
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname manhattanplot
#' @aliases manhattanplot,compareClusterResult,ANY-method
#' @exportMethod manhattanplot
setMethod(
    "manhattanplot",
    signature(x = "compareClusterResult"),
    function(
        x,
        color = "p.adjust",
        showCategory = 5,
        split = NULL,
        font.size = 12,
        title = "",
        size = "Count",
        includeAll = TRUE,
        label_format = 30,
        ...
    ) {
        manhattanplot.compareClusterResult(
            x,
            colorBy = color,
            showCategory = showCategory,
            size = size,
            includeAll = includeAll,
            split = split,
            font.size = font.size,
            title = title,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname manhattanplot
#' @exportMethod manhattanplot
#' @aliases manhattanplot,enrichResultList,ANY-method
setMethod(
    "manhattanplot",
    signature(x = "enrichResultList"),
    function(
        x,
        color = "p.adjust",
        showCategory = 5,
        size = "Count",
        split = NULL,
        font.size = 12,
        title = "",
        label_format = 30,
        ...
    ) {
        manhattanplot.enrichResult(
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname manhattanplot
#' @exportMethod manhattanplot
#' @aliases manhattanplot,gseaResultList,ANY-method
setMethod(
    "manhattanplot",
    signature(x = "gseaResultList"),
    function(
        x,
        color = "p.adjust",
        showCategory = 5,
        size = "Count",
        split = NULL,
        font.size = 12,
        title = "",
        label_format = 30,
        ...
    ) {
        manhattanplot.enrichResult(
            x = x,
            color = color,
            showCategory = showCategory,
            size = size,
            split = split,
            font.size = font.size,
            title = title,
            label_format = label_format,
            ...
        )
    }
)

#' @rdname manhattanplot
#' @exportMethod manhattanplot
#' @aliases manhattanplot,list,ANY-method
setMethod(
    "manhattanplot",
    signature(x = "list"),
    function(
        x,
        color = "p.adjust",
        showCategory = 5,
        size = "Count",
        split = NULL,
        font.size = 12,
        title = "",
        label_format = 30,
        ...
    ) {
        if (all(sapply(x, function(i) inherits(i, "enrichResult") || inherits(i, "gseaResult")))) {
            class(x) <- "enrichResultList"
            manhattanplot(
                x = x,
                color = color,
                showCategory = showCategory,
                size = size,
                split = split,
                font.size = font.size,
                title = title,
                label_format = label_format,
                ...
            )
        } else {
            stop("all elements in the list should be enrichResult or gseaResult objects")
        }
    }
)

#' Internal helper function for manhattan build
#' @noRd

.get_ontology <- function(x) {
    if ("ontology" %in% methods::slotNames(x) && length(x@ontology) > 0 && x@ontology != "") {
        return(x@ontology)
    }
    if ("fun" %in% methods::slotNames(x) && length(x@fun) > 0 && x@fun != "") {
        if (x@fun == "enrichGO" || x@fun == "gseGO") {
           return("GO")
        }
        res <- gsub("enrich", "", x@fun) 
        res <- gsub("gse", "", res)
        return(res)
    }
    return("Enrichment")
}

.prep_manhattan_df <- function(df, colorBy) {
    if (nrow(df) == 0) return(list(df = df))
    
    grp_col <- "ONTOLOGY"
    if (!"ONTOLOGY" %in% colnames(df)) {
        if ("Category" %in% colnames(df)) {
            grp_col <- "Category"
        } else {
            df$ONTOLOGY <- "Enrichment"
        }
    }
    
    unique_terms <- unique(df[, c("ID", grp_col)])
    unique_terms <- unique_terms[order(unique_terms[[grp_col]], unique_terms$ID), ]
    
    unique_terms$x_pos <- NA
    grps <- unique(unique_terms[[grp_col]])
    
    current_x <- 0
    ticks <- numeric(length(grps))
    gap <- max(1, nrow(unique_terms) * 0.05)
    
    for (i in seq_along(grps)) {
        grp <- grps[i]
        idx <- which(unique_terms[[grp_col]] == grp)
        n <- length(idx)
        unique_terms$x_pos[idx] <- current_x + (1:n)
        ticks[i] <- current_x + n / 2
        current_x <- current_x + n + gap
    }
    
    df <- merge(df, unique_terms, by = c("ID", grp_col))
    
    # Calculate y
    df$y <- -log10(df[[colorBy]])
    
    list(df = df, ticks = ticks, grps = grps, grp_col = grp_col)
}

.manhattanplot_internal <- function(df, ticks, grps, grp_col, hl_df, size, colorBy, label_func, font.size, title, size_range = c(3, 8)) {
    p <- ggplot(df, aes(x = .data$x_pos, y = .data$y)) +
        geom_point(aes(size = .data[[size]], fill = .data[[grp_col]]), shape = 21, alpha = 0.8) +
        scale_x_continuous(breaks = ticks, labels = grps) +
        ylab(paste0("-log10(", colorBy, ")")) +
        xlab(NULL) +
        ggtitle(title) +
        theme_dose(font.size) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    
    if (size == "Count" && !is.null(df[[size]])) {
        tryCatch({
            size_break <- pretty(df[[size]], n = 4)
            p <- p + scale_size(range = size_range, breaks = size_break)
        }, error = function(e) {
            p <- p + scale_size(range = size_range)
        })
    } else {
        p <- p + scale_size(range = size_range)
    }
    
    if (nrow(hl_df) > 0) {
        require_suggested('ggrepel', 'for labeling in `manhattanplot()`.')
        hl_df$label <- label_func(hl_df$Description)
        p <- p + ggrepel::geom_text_repel(
            data = hl_df,
            aes(x = .data$x_pos, y = .data$y, label = .data$label),
            size = font.size / 3,
            min.segment.length = 0,
            box.padding = 0.5,
            show.legend = FALSE
        )
    }
    
    class(p) <- c("enrichplotManhattan", class(p))
    return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous xlab ylab ggtitle theme element_text element_blank scale_size theme_void
#' @importFrom utils head
manhattanplot.enrichResult <- function(
    x,
    color = "p.adjust",
    showCategory = 5,
    size = "Count",
    split = NULL,
    font.size = 12,
    title = "",
    label_format = 30,
    ...
) {
    colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
    size <- normalize_measure_var(size)

    if (inherits(x, c("enrichResultList", "gseaResultList"))) {
        ldf <- lapply(x, as.data.frame)
        n_all <- sum(sapply(ldf, nrow))
        if (n_all == 0) return(ggplot() + theme_void())
        
        ldf <- lapply(seq_along(x), function(i) {
            obj <- x[[i]]
            df_i <- fortify(obj, showCategory = nrow(as.data.frame(obj)))
            if (!"ONTOLOGY" %in% colnames(df_i)) {
               df_i$ONTOLOGY <- .get_ontology(obj)
            }
            return(df_i)
        })
        names(ldf) <- names(x)
        df <- dplyr::bind_rows(ldf, .id = "category")
        df$category <- factor(df$category, levels = names(x))
    } else {
        n_all <- nrow(as.data.frame(x))
        if (n_all == 0) return(ggplot() + theme_void())
        df <- fortify(x, showCategory = n_all, split = split)
        if (!"ONTOLOGY" %in% colnames(df)) {
           df$ONTOLOGY <- .get_ontology(x)
        }
    }

    res <- .prep_manhattan_df(df, colorBy)
    df <- res$df
    
    if (!is.null(showCategory) && showCategory > 0) {
        df_ord <- df[order(df$y, decreasing = TRUE), ]
        hl_df <- head(df_ord, showCategory)
    } else {
        hl_df <- df[0, ]
    }

    label_func <- .label_format(label_format)

    p <- .manhattanplot_internal(
        df = df,
        ticks = res$ticks,
        grps = res$grps,
        grp_col = res$grp_col,
        hl_df = hl_df,
        size = size,
        colorBy = colorBy,
        label_func = label_func,
        font.size = font.size,
        title = title
    )
    
    return(p)
}

#' @importFrom ggplot2 facet_grid
manhattanplot.compareClusterResult <- function(
    x,
    colorBy = "p.adjust",
    showCategory = 5,
    size = "Count",
    split = NULL,
    includeAll = TRUE,
    font.size = 12,
    title = "",
    label_format = 30,
    facet = "Cluster",
    strip_width = 15,
    ...
) {
    size <- normalize_measure_var(size)

    if (!is.null(facet) && facet == "intersect") {
        x <- append_intersect(x)
    }

    n_all <- nrow(as.data.frame(x))
    if (n_all == 0) return(ggplot() + theme_void())
    
    df <- fortify(x, showCategory = n_all, includeAll = includeAll, split = split)
    if (!"ONTOLOGY" %in% colnames(df)) {
       df$ONTOLOGY <- .get_ontology(x)
    }

    # In single enrich we didn't use `colorBy`, we used `color`. For compare we mapped it to colorBy.
    colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue"))

    res <- .prep_manhattan_df(df, colorBy)
    df <- res$df

    if (!is.null(showCategory) && showCategory > 0) {
        df_ord <- df[order(df$y, decreasing = TRUE), ]
        hl_df <- do.call(rbind, by(df_ord, df_ord[[facet]], head, n = showCategory))
    } else {
        hl_df <- df[0, ]
    }

    label_func <- .label_format(label_format)

    p <- .manhattanplot_internal(
        df = df,
        ticks = res$ticks,
        grps = res$grps,
        grp_col = res$grp_col,
        hl_df = hl_df,
        size = size,
        colorBy = colorBy,
        label_func = label_func,
        font.size = font.size,
        title = title
    )
    
    if (!is.null(facet)) {
        p <- p +
            facet_grid(
                stats::reformulate(".", response = facet),
                scales = "free_y",
                space = 'fixed',
                switch = 'y',
                labeller = ggplot2::label_wrap_gen(strip_width)
            ) +
            theme(strip.text = element_text(size = 14))
    }
    
    return(p)
}
