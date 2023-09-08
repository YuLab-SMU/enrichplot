##' @rdname ssplot
##' @exportMethod ssplot
setMethod("ssplot", signature(x = "enrichResult"),
          function(x, showCategory = 30,  ...) {
              ssplot.enrichResult(x, showCategory = showCategory, ...)
          })

##' @rdname ssplot
##' @exportMethod ssplot
setMethod("ssplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, ...) {
              ssplot.enrichResult(x, showCategory = showCategory, ...)
          })

##' @rdname ssplot
##' @exportMethod ssplot
setMethod("ssplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 30, ...) {
              ssplot.compareClusterResult(x, showCategory = showCategory,
                                         ...)
          })




##' @rdname ssplot
##' @param drfun The function used for dimension reduction,
##' e.g. stats::cmdscale (the default), vegan::metaMDS, or ape::pcoa.
##' @param dr.params list, the parameters of tidydr::dr. one of 'category', 'group', 'all' and 'none'.
##' @param ... additional parameters
##' 
##' additional parameters can refer the following parameters.
##'     \itemize{
##'        \item \code{coords} a data.frame with two columns: 'x' for X-axis coordinate and 'y' for Y-axis coordinate.
##'        \item \code{color} Variable that used to color enriched terms, e.g. 'pvalue','p.adjust' or 'qvalue'.
##'         the starting position of each text label. 
##'        \item \code{cex_line} Scale of line width.
##'        \item \code{min_edge} The minimum similarity threshold for whether 
##'         two nodes are connected, should between 0 and 1, default value is 0.2.
##'        \item \code{cex_label_category} Scale of category node label size.
##'        \item \code{cex_category} Number indicating the amount by which plotting category
##'         nodes should be scaled relative to the default.
##'        \item \code{shadowtext} a logical value, whether to use shadow font. 

##'        \item \code{label_style} style of group label, one of "shadowtext" and "ggforce".
##'        \item \code{repel whether} to correct the position of the label. Defaults to FALSE.
##'        \item \code{group_legend} Logical, if TRUE, the grouping legend will be displayed.
##'         The default is FALSE.
##'        \item \code{cex_label_group} Numeric, scale of group labels size, the default value is 1.
##'        \item \code{nWords} Numeric, the number of words in the cluster tags, the default value is 4.
##'        \item \code{label_format} a numeric value sets wrap length, alternatively a
##'         custom function to format axis labels.
##'        \item \code{clusterFunction} function of Clustering method, such as stats::kmeans(the default),
##'         cluster::clara, cluster::fanny or cluster::pam.
##'        \item \code{nCluster} Numeric, the number of clusters, 
##'         the default value is square root of the number of nodes.
##'     }
##'
##' additional parameters can refer the emapplot function: \link{emapplot}.
ssplot.enrichResult <- function(x, showCategory = 30,
                                drfun = NULL,
                                with_edge = FALSE,
                                dr.params = list(),
                                group_category = TRUE,
                                node_label  = "group",
                                ...) {
    if (is.null(drfun)) {
        drfun = stats::cmdscale
        dr.params = list(eig = TRUE)
    }
    if (is.character(drfun)) {
        drfun <- eval(parse(text=drfun))
    }

    drResult <- get_drResult(x = x, showCategory = showCategory, 
        drfun = drfun, dr.params = dr.params)
    coords <- drResult$drdata[, c(1, 2)]
    colnames(coords) <- c("x", "y")
    rownames(coords) <- attr(drResult$data, "Labels")
    p <- emapplot(x = x, showCategory = showCategory,
                  coords = coords,
                  with_edge = with_edge,
                  group_category = group_category,
                  node_label = node_label,
                  ...)

    ## Set axis label according to drfun
    p <- adj_axis(p = p, drResult = drResult)
    return(p + theme_classic())
}




##' @rdname ssplot
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 coord_equal
##' @importClassesFrom DOSE compareClusterResult
##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis,
##' the default value is 0.0125.
##' @importFrom stats setNames
ssplot.compareClusterResult <- function(x, showCategory = 30,
                                        split = NULL, pie = "equal",
                                        drfun = NULL,
                                        with_edge = FALSE,
                                        cex_pie2axis = 0.0125, 
                                        dr.params = list(), 
                                        group_category = TRUE,
                                        node_label  = "group", ...) {
    if (is.null(drfun)) {
        drfun = stats::cmdscale
        dr.params = list(eig = TRUE)
    }

    if (is.character(drfun)) {
        drfun <- eval(parse(text=drfun))
    }

    drResult <- get_drResult(x = x, showCategory = showCategory, 
        split = split, pie = pie, drfun = drfun, dr.params = dr.params)
    coords <- drResult$drdata[, c(1, 2)]
    colnames(coords) <- c("x", "y")
    rownames(coords) <- attr(drResult$data, "Labels")
    p <- emapplot(x, showCategory = showCategory,
                  coords = coords,
                  split = split, pie = pie,
                  with_edge = with_edge,
                  cex_pie2axis = cex_pie2axis, 
                  group_category = group_category,
                  node_label = node_label, ...)
    ## Set axis label according to the method parameter
    p <- adj_axis(p = p, drResult = drResult)
    return(p + theme_classic())
}


##' Get a distance matrix
##'
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display.
##' @param split separate result by 'category' variable.
##' @param pie proportion of clusters in the pie chart.
##' @noRd
build_dist <- function(x, showCategory, split = NULL, pie = NULL) {
    sim = get_pairwise_sim(x = x, showCategory = showCategory,
        split = split, pie = pie)
    ## If the similarity between the two terms is 1,
    ## an error will be reported in some method, so fine-tuning.
    sim[which(sim == 1)] <- 0.99999
    for (i in seq_len(nrow(sim))) sim[i, i] <- 1 
    stats::as.dist(1- sim)
}


##' Get a similarity matrix
##'
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display.
##' @param split separate result by 'category' variable.
##' @param pie proportion of clusters in the pie chart.
##' @noRd
get_pairwise_sim <- function(x, showCategory, split = NULL, pie = NULL) {
    # if (class(x) == "compareClusterResult") {
    if (inherits(x, "compareClusterResult")) {
        # y <- get_selected_category(showCategory, enrichResult, split)
        y <- fortify(model = x, showCategory = showCategory,
                     includeAll = TRUE, split = split)
        y$Cluster <- sub("\n.*", "", y$Cluster)
        keep <- rownames(prepare_pie_category(y, pie=pie))
    } else {
        n <- update_n(x, showCategory)
        if (is.numeric(n)) {
            keep <- seq_len(n)
        } else {
            keep <- match(n, rownames(x@termsim))
        }
    }
    if (length(keep) == 0) {
        stop("no enriched term found...")
    }
    fill_termsim(x, keep)
}



##' Adjust axis label according to the dimension reduction method
##'
##' @param p ggplot2 object
##' @param drs dimension reduction result
##' @noRd
adj_axis <- function(p, drResult) {
    title = NULL
    eigenvalue <- drResult$eigenvalue
    if (!is.null(eigenvalue)) {
        xlab = paste("Dimension1 (", format(100 * eigenvalue[1] / sum(eigenvalue), digits=4), "%)", sep = "")
        ylab = paste("Dimension2 (", format(100 * eigenvalue[2] / sum(eigenvalue), digits=4), "%)", sep = "")
    } else {
        xlab = "Dimension1"
        ylab = "Dimension2"
        if (!is.null(drResult$stress)) {
            title = paste("stress = ", drResult$stress, sep = "")
        }
    }
    p <- p +  labs(x = xlab, y = ylab, title = title)
    return(p)
}

##' Get the result of dimension reduction
##'
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display.
##' @param split separate result by 'category' variable.
##' @param pie proportion of clusters in the pie chart.
##' @param drfun The function used for dimension reduction.
##' @param dr.params list, the parameters of tidydr::dr.
##' @importFrom rlang check_installed
##' @noRd
get_drResult <- function(x, showCategory, split = NULL, pie = NULL, drfun, dr.params) {
    distance_mat <- build_dist(x = x, showCategory = showCategory, split = split, pie = pie)
    check_installed('tidydr', 'for `get_drResult()`, which is an internal function.')
    drResult <- do.call(tidydr::dr, c(list(data = distance_mat, fun = drfun), dr.params))
    wrongMessage <- paste("Wrong drfun parameter or unsupported",
         "dimensionality reduction method;",
         "set to default `drfun = 'stats::cmdscale'`")
    if (is.null(drResult$drdata)) {
        message(wrongMessage)
        drResult <- tidydr::dr(distance_mat, stats::cmdscale, eig = TRUE)
    }
    drResult
}
