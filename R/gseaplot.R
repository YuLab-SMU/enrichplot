#' @rdname gseaplot
#' @exportMethod gseaplot
setMethod(
    "gseaplot",
    signature(x = "gseaResult"),
    function(
        x,
        geneSetID,
        by = "all",
        title = "",
        color = 'black',
        color.line = "green",
        color.vline = "#FA5860",
        ...
    ) {
        gseaplot.gseaResult(
            x,
            geneSetID = geneSetID,
            by = by,
            title = title,
            color = color,
            color.line = color.line,
            color.vline = color.vline,
            ...
        )
    }
)

#' @rdname gseaplot
#' @param layer Optional `mnsea` layer. When `NULL`, use collapsed scores.
#' @exportMethod gseaplot
setMethod(
    "gseaplot",
    signature(x = "mnseaResult"),
    function(
        x,
        geneSetID,
        by = "all",
        title = "",
        color = "black",
        color.line = "green",
        color.vline = "#FA5860",
        layer = NULL,
        ...
    ) {
        gseaplot.mnseaResult(
            x,
            geneSetID = geneSetID,
            by = by,
            title = title,
            color = color,
            color.line = color.line,
            color.vline = color.vline,
            layer = layer,
            ...
        )
    }
)

#' @rdname gseaplot
#' @param color color of line segments
#' @param color.line color of running enrichment score line
#' @param color.vline color of vertical line indicating the
#' maximum/minimal running enrichment score
#' @return ggplot2 object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplotGrob
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 rel
#' @importFrom aplot plot_list
#' @author Guangchuang Yu
gseaplot.gseaResult <- function(
    x,
    geneSetID,
    by = "all",
    title = "",
    color = 'black',
    color.line = "green",
    color.vline = "#FA5860",
    ...
) {
    by <- match.arg(by, c("runningScore", "preranked", "all"))
    gsdata <- gsInfo(x, geneSetID)
    p <- ggplot(gsdata, aes(x = .data$x)) +
        theme_dose() +
        xlab("Position in the Ranked List of Genes")
    if (by == "runningScore" || by == "all") {
        p.res <- p +
            geom_linerange(
                aes(ymin = .data$ymin, ymax = .data$ymax),
                color = color
            )
        p.res <- p.res +
            geom_line(
                aes(y = .data$runningScore),
                color = color.line,
                linewidth = 1
            )
        enrichmentScore <- x@result[geneSetID, "enrichmentScore"]
        es.df <- data.frame(
            es = which.min(abs(p$data$runningScore - enrichmentScore))
        )
        p.res <- p.res +
            geom_vline(
                data = es.df,
                aes(xintercept = .data$es),
                colour = color.vline,
                linetype = "dashed"
            )
        p.res <- p.res + ylab("Running Enrichment Score")
        p.res <- p.res + geom_hline(yintercept = 0)
    }
    if (by == "preranked" || by == "all") {
        df2 <- data.frame(x = which(p$data$position == 1))
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p +
            geom_segment(
                data = df2,
                aes(x = .data$x, xend = .data$x, y = .data$y, yend = 0),
                color = color
            )
        p.pos <- p.pos +
            ylab("Ranked List Metric") +
            xlim(0, length(p$data$geneList))
    }
    if (by == "runningScore") {
        return(p.res + ggtitle(title))
    }
    if (by == "preranked") {
        return(p.pos + ggtitle(title))
    }

    p.pos <- p.pos +
        xlab(NULL) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p.pos <- p.pos +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
    #plot_list(gglist =  list(p.pos, p.res), ncol=1)

    aplot::gglist(gglist = list(p.pos, p.res), ncol = 1)
}

gseaplot.mnseaResult <- function(
    x,
    geneSetID,
    by = "all",
    title = "",
    color = "black",
    color.line = "green",
    color.vline = "#FA5860",
    layer = NULL,
    ...
) {
    by <- match.arg(by, c("runningScore", "preranked", "all"))
    gsdata <- gsInfo(x, geneSetID = geneSetID, layer = layer)
    p <- ggplot(gsdata, aes(x = .data$x)) +
        theme_dose() +
        xlab("Position in the Ranked List of Genes")
    if (by == "runningScore" || by == "all") {
        p.res <- p +
            geom_linerange(
                aes(ymin = .data$ymin, ymax = .data$ymax),
                color = color
            )
        p.res <- p.res +
            geom_line(
                aes(y = .data$runningScore),
                color = color.line,
                linewidth = 1
            )
        enrichmentScore <- x@result[gsdata$ID[1], "enrichmentScore"]
        es.df <- data.frame(
            es = which.min(abs(p$data$runningScore - enrichmentScore))
        )
        p.res <- p.res +
            geom_vline(
                data = es.df,
                aes(xintercept = .data$es),
                colour = color.vline,
                linetype = "dashed"
            )
        p.res <- p.res + ylab("Running Enrichment Score")
        p.res <- p.res + geom_hline(yintercept = 0)
    }
    if (by == "preranked" || by == "all") {
        df2 <- data.frame(x = which(p$data$position == 1))
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p +
            geom_segment(
                data = df2,
                aes(x = .data$x, xend = .data$x, y = .data$y, yend = 0),
                color = color
            )
        p.pos <- p.pos +
            ylab("Ranked List Metric") +
            xlim(0, length(p$data$geneList))
    }
    if (by == "runningScore") {
        return(p.res + ggtitle(title))
    }
    if (by == "preranked") {
        return(p.pos + ggtitle(title))
    }

    p.pos <- p.pos +
        xlab(NULL) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p.pos <- p.pos +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5, size = rel(2)))

    aplot::gglist(gglist = list(p.pos, p.res), ncol = 1)
}

get_mnsea_ranked_scores <- function(object, layer = NULL) {
    if (is.null(layer)) {
        geneList <- object@collapsed_scores
    } else {
        layer <- unique(as.character(layer))
        if (length(layer) != 1) {
            yulab.utils::yulab_abort(
                "`layer` must be `NULL` or a single layer for `gseaplot()`."
            )
        }
        keep_layer <- resolve_mnsea_layers(names(object@layer_scores), layer = layer)
        if (length(keep_layer) == 0) {
            yulab.utils::yulab_abort(
                "No mnsea scores available for the selected `layer`."
            )
        }
        geneList <- object@layer_scores[[keep_layer[1]]]
    }

    geneList <- geneList[!is.na(geneList)]
    if (length(geneList) == 0 || is.null(names(geneList))) {
        yulab.utils::yulab_abort("No ranked mnsea scores available for plotting.")
    }
    geneList
}

resolve_mnsea_gseaplot_id <- function(object, geneSetID) {
    if (is.null(geneSetID)) {
        return(default_mnsea_pathway_id(object))
    }

    geneSetID <- unname(geneSetID)
    if (length(geneSetID) != 1) {
        yulab.utils::yulab_abort(
            "`gseaplot.mnseaResult()` only supports a single pathway."
        )
    }

    result_df <- .result_data(object)
    if (is.numeric(geneSetID)) {
        idx <- as.integer(geneSetID)
        if (is.na(idx) || idx < 1 || idx > nrow(result_df)) {
            yulab.utils::yulab_abort("`geneSetID` is out of bounds.")
        }
        return(as.character(result_df$ID[idx]))
    }

    idx <- resolve_term_rows(object, geneSetID)
    if (length(idx) == 0) {
        yulab.utils::yulab_abort(
            "No mnsea pathway matched the supplied `geneSetID`."
        )
    }
    as.character(result_df$ID[idx[1]])
}

get_mnsea_geneSet <- function(object, geneSetID, layer = NULL) {
    feature_df <- fortify_mnsea_contribution(
        object,
        level = "feature",
        pathway_id = geneSetID,
        layer = layer
    )
    geneSet <- unique(as.character(feature_df$Feature))
    geneSet <- geneSet[!is.na(geneSet) & nzchar(geneSet)]

    if (length(geneSet) == 0 && geneSetID %in% names(object@geneSets)) {
        geneSet <- unique(as.character(object@geneSets[[geneSetID]]))
        geneSet <- geneSet[!is.na(geneSet) & nzchar(geneSet)]
    }

    if (length(geneSet) == 0) {
        yulab.utils::yulab_abort(
            "No mnsea features available for the selected pathway/layer."
        )
    }
    geneSet
}

#' @method gsInfo mnseaResult
#' @export
gsInfo.mnseaResult <- function(object, geneSetID, layer = NULL, ...) {
    geneSetID <- resolve_mnsea_gseaplot_id(object, geneSetID)
    geneList <- get_mnsea_ranked_scores(object, layer = layer)
    geneSet <- get_mnsea_geneSet(object, geneSetID = geneSetID, layer = layer)
    geneSet <- intersect(geneSet, names(geneList))
    if (length(geneSet) == 0) {
        yulab.utils::yulab_abort(
            "No overlap between mnsea features and ranked scores for the selected pathway/layer."
        )
    }

    exponent <- object@params[["exponent"]]
    if (length(exponent) == 0 || is.na(exponent)) {
        exponent <- 1
    }

    df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore)) / 20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
    } else {
        df$gene <- object@gene2Symbol[names(geneList)]
    }

    description <- as.character(get_term_labels(object, geneSetID)[[1]])
    df$ID <- geneSetID
    df$Description <- description
    df$layer <- if (is.null(layer)) "collapsed" else as.character(layer)
    df
}


#' extract gsea result of selected geneSet
#'
#'
#' @title gsInfo
#' @param object gseaResult object
#' @param geneSetID gene set ID
#' @param ... additional arguments passed to methods (e.g. `layer` for `mnseaResult`)
#' @return data.frame
#' @author Guangchuang Yu
gsInfo <- function(object, geneSetID, ...) {
    UseMethod("gsInfo")
}

#' @method gsInfo gseaResult
#' @export
gsInfo.gseaResult <- function(object, geneSetID, ...) {
    geneList <- object@geneList

    if (is.numeric(geneSetID)) {
        geneSetID <- object@result[geneSetID, "ID"]
    }

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    if (length(exponent) == 0 || is.na(exponent)) {
        exponent <- 1
    }
    df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore)) / 20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
    } else {
        df$gene <- object@gene2Symbol[names(geneList)]
    }

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}


get_gsdata <- function(x, geneSetID, ...) {
    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID, ...)
        return(gsdata)
    }

    lapply(geneSetID, gsInfo, object = x, ...) |>
        yulab.utils::rbindlist()
}

#' Horizontal plot for GSEA result
#'
#'
#' @title hplot
#' @param x gseaResult object
#' @param geneSetID gene set ID
#' @param layer Optional `mnsea` layer. When `NULL`, use collapsed scores.
#' @return horizontal plot
#' @export
#' @author Guangchuang Yu
hplot <- function(x, geneSetID, layer = NULL) {
    if (!inherits(x, "gseaResult")) {
        stop("hplot only work for GSEA result")
    }

    gsdata <- get_gsdata(x, geneSetID, layer = layer)

    ggplot(gsdata, aes(.data$x, .data$runningScore)) +
        ggplot2::geom_ribbon(
            aes(
                ymin = pmin(0, .data$runningScore),
                ymax = pmax(0, .data$runningScore),
                fill = ifelse(.data$runningScore >= 0, "up", "down")
            ),
            alpha = 0.5
        ) +
        geom_line(color = "grey30", linewidth = 0.3) +
        geom_hline(yintercept = 0, color = "grey50", linewidth = 0.2) +
        facet_grid(Description ~ .) +
        scale_fill_manual(
            values = c(up = "#FA5860", down = "#2166AC"),
            name = NULL
        ) +
        theme_minimal() +
        ggfun::theme_noyaxis() +
        theme(
            panel.spacing.y = unit(0, "lines"),
            strip.text.y = element_text(angle = 0),
            legend.position = "none",
            panel.border = element_blank(),
            panel.grid = element_blank(),
        ) +
        xlab(NULL) +
        ylab(NULL)
}

#' GSEA plot that mimic the plot generated by broad institute's GSEA software
#'
#'
#' @title gseaplot2
#' @param x gseaResult object
#' @param geneSetID gene set ID
#' @param title plot title
#' @param color color of running enrichment score line
#' @param base_size base font size
#' @param rel_heights relative heights of subplots
#' @param subplots which subplots to be displayed
#' @param pvalue_table whether add pvalue table
#' @param pvalue_table_columns selected columns to be plotted in the `pvalue_table`
#' @param pvalue_table_rownames selected column as the rownames of the `pvalue_table`. If set to NULL, no rownames will be displayed.
#' @param ES_geom geom for plotting running enrichment score,
#' one of 'line' or 'dot'
#' @param layer Optional `mnsea` layer. When `NULL`, use collapsed scores.
#' @return plot
#' @export
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 annotation_custom
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @author Guangchuang Yu
gseaplot2 <- function(
    x,
    geneSetID,
    title = "",
    color = "green",
    base_size = 11,
    rel_heights = c(1.5, .5, 1),
    subplots = 1:3,
    pvalue_table = FALSE,
    pvalue_table_columns = c("pvalue", "p.adjust"),
    pvalue_table_rownames = "Description",
    ES_geom = "line",
    layer = NULL
) {
    ES_geom <- match.arg(ES_geom, c("line", "dot"))

    geneList <- position <- NULL ## to satisfy codetool

    gsdata <- get_gsdata(x, geneSetID, layer = layer)

    p <- ggplot(gsdata, aes(x = .data$x)) +
        xlab(NULL) +
        theme_classic(base_size) +
        theme(
            panel.grid.major = element_line(colour = "grey92"),
            panel.grid.minor = element_line(colour = "grey92"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
        ) +
        scale_x_continuous(expand = c(0, 0))

    if (ES_geom == "line") {
        es_layer <- geom_line(
            aes(y = .data$runningScore, color = .data$Description),
            linewidth = 1
        )
    } else {
        es_layer <- geom_point(
            aes(y = .data$runningScore, color = .data$Description),
            size = 1,
            data = subset(gsdata, position == 1)
        )
    }

    p.res <- p +
        es_layer +
        theme(
            legend.position = "inside",
            legend.position.inside = c(.8, .8),
            legend.title = element_blank(),
            legend.background = element_rect(fill = "transparent")
        )

    p.res <- p.res +
        ylab("Running Enrichment Score") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
        )

    # Vectorized ymin/ymax assignment
    terms <- unique(gsdata$Description)
    term_indices <- match(gsdata$Description, terms) - 1
    idx <- which(gsdata$ymin != 0)
    gsdata[idx, "ymin"] <- term_indices[idx]
    gsdata[idx, "ymax"] <- term_indices[idx] + 1
    p2 <- ggplot(gsdata, aes(x = .data$x)) +
        geom_linerange(aes(
            ymin = .data$ymin,
            ymax = .data$ymax,
            color = .data$Description
        )) +
        xlab(NULL) +
        ylab(NULL) +
        theme_classic(base_size) +
        theme(
            legend.position = "none",
            plot.margin = margin(t = -.1, b = 0, unit = "cm"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_blank()
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))

    if (length(geneSetID) == 1) {
        ## geneList <- gsdata$geneList
        ## j <- which.min(abs(geneList))
        ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
        ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]

        ## v <- sort(c(v1, v2))
        ## inv <- findInterval(geneList, v)

        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) {
            inv <- inv + 1
        }

        col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * .3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(
            ymin = ymin,
            ymax = yy,
            xmin = xmin,
            xmax = xmax,
            col = col[unique(inv)]
        )
        p2 <- p2 +
            geom_rect(
                aes(
                    xmin = .data$xmin,
                    xmax = .data$xmax,
                    ymin = .data$ymin,
                    ymax = .data$ymax,
                    fill = I(col)
                ),
                data = d,
                alpha = .9,
                inherit.aes = FALSE
            )
    }

    ## p2 <- p2 +
    ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
    ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
    ## theme(legend.position="none") +
    ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))

    df2 <- p$data #data.frame(x = which(p$data$position == 1))
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p +
        geom_segment(
            data = df2,
            aes(x = .data$x, xend = .data$x, y = .data$y, yend = 0),
            color = "grey"
        )
    p.pos <- p.pos +
        ylab("Ranked List Metric") +
        xlab("Rank in Ordered Dataset") +
        theme(
            plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm")
        )

    if (!is.null(title) && !is.na(title) && title != "") {
        p.res <- p.res + ggtitle(title)
    }

    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        } else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }

    if (pvalue_table) {
        pd <- x[geneSetID, pvalue_table_columns]
        # pd <- pd[order(pd[,1], decreasing=FALSE),]
        if (is.null(pvalue_table_rownames)) {
            rows <- NULL
        } else {
            # rownames(pd) <- pd$Description
            if (length(pvalue_table_rownames) != 1) {
                stop(
                    "the length of `pvalue_table_rownames` should be equal to 1"
                )
            }

            rows <- x[geneSetID, pvalue_table_rownames]
        }

        # pd <- round(pd, 4)
        for (i in seq_len(ncol(pd))) {
            pd[, i] <- format(pd[, i], digits = 4)
        }
        tp <- tableGrob2(d = pd, p = p.res, rows = rows)

        p.res <- p.res +
            theme(legend.position = "none") +
            annotation_custom(
                tp,
                xmin = quantile(p.res$data$x, .5),
                xmax = quantile(p.res$data$x, .95),
                ymin = quantile(p.res$data$runningScore, .75),
                ymax = quantile(p.res$data$runningScore, .9)
            )
    }

    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] +
        theme(
            axis.line.x = element_line(),
            axis.ticks.x = element_line(),
            axis.text.x = element_text()
        )

    if (length(subplots) == 1) {
        return(
            plotlist[[1]] +
                theme(
                    plot.margin = margin(
                        t = .2,
                        r = .2,
                        b = .2,
                        l = .2,
                        unit = "cm"
                    )
                )
        )
    }

    if (length(rel_heights) > length(subplots)) {
        rel_heights <- rel_heights[subplots]
    }

    # aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
    aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights)
}


#' plot ranked list of genes with running enrichment score as bar height
#'
#'
#' @title gsearank
#' @param x gseaResult object
#' @param geneSetID gene set ID
#' @param title plot title
#' @param output one of 'plot' or 'table' (for exporting data)
#' @param layer Optional `mnsea` layer. When `NULL`, use collapsed scores.
#' @return ggplot object
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 theme_minimal
#' @export
#' @author Guangchuang Yu
gsearank <- function(x, geneSetID, title = "", output = "plot", layer = NULL) {
    output <- match.arg(output, c("plot", "table"))

    position <- NULL
    gsdata <- gsInfo(x, geneSetID, layer = layer)
    gsdata <- subset(gsdata, position == 1)

    if (output == "table") {
        res <- gsdata[, c("gene", "x", "runningScore")]
        if (x[geneSetID, "NES"] > 0) {
            res$core <- "NO"
            res$core[1:which.max(gsdata$runningScore)] <- "YES"
        } else {
            res$core <- "NO"
            res$core[which.min(gsdata$runningScore):nrow(res)] <- "YES"
        }
        names(res) <- c(
            "gene",
            "rank in geneList",
            "running ES",
            "core enrichment"
        )
        rownames(res) <- NULL
        return(res)
    }

    p <- ggplot(gsdata, aes(x = .data$x, y = .data$runningScore)) +
        geom_segment(aes(xend = .data$x, yend = 0)) +
        ggtitle(title) +
        xlab("Position in the Ranked List of Genes") +
        ylab("Running Enrichment Score") +
        theme_minimal()
    return(p)
}


#' label genes in running score plot
#'
#'
#' @title geom_gsea_gene
#' @param genes selected genes to be labeled
#' @param mapping aesthetic mapping, default is NULL
#' @param geom geometric layer to plot the gene labels, default is geom_text
#' @param ... additional parameters passed to the 'geom'
#' @param geneSet choose which gene set(s) to be label if the plot contains multiple gene sets
#' @return ggplot object
#' @importFrom rlang .data
#' @export
#' @author Guangchuang Yu
geom_gsea_gene <- function(
    genes,
    mapping = NULL,
    geom = ggplot2::geom_text,
    ...,
    geneSet = NULL
) {
    default_mapping <- aes(
        x = .data$x,
        y = .data$runningScore,
        label = .data$gene
    )
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }
    if (is.null(geneSet)) {
        data <- ggtree::td_filter(.data$gene %in% genes)
    } else {
        data <- ggtree::td_filter(
            .data$gene %in% genes & .data$Description %in% geneSet
        )
    }

    geom(mapping = mapping, data = data, ...)
}
