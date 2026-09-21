#' Get the similarity matrix
#'
#' @param y A data.frame of enrichment result
#' @param geneSets A list, the names of geneSets are term ids,
#' and every element is a vector of genes.
#' @param method Method of calculating the similarity between nodes,
#' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
#' "JC" (Jaccard similarity coefficient) methods
#' @param semData GOSemSimDATA object
#' @noRd
get_similarity_matrix <- function(y, geneSets, method, semData = NULL) {
    id <- y[, "ID"]
    geneSets <- geneSets[id]
    y_id <- unlist(strsplit(y$ID[1], ":"))[1]
    ## Choose the method to calculate the similarity
    if (method == "JC") {
        w <- .cal_jc_similarity(geneSets, id = id, name = y$Description)
        rownames(y) <- y$ID
        rownames(w) <- colnames(w) <- unname(get_term_labels(y, id))
        return(w)
    }

    if (y_id == "GO") {
        if (is.null(semData)) {
            stop(
                "The semData parameter is missing,
                and it can be obtained through godata function in GOSemSim package."
            )
        }
        w <- GOSemSim::mgoSim(
            id,
            id,
            semData = semData,
            measure = method,
            combine = NULL
        )
    }

    if (y_id == "DOID") {
        w <- DOSE::doSim(id, id, measure = method)
    }
    rownames(y) <- y$ID
    rownames(w) <- colnames(w) <- unname(get_term_labels(y, colnames(w)))
    return(w)
}


#' Check whether the similarity matrix exists
#'
#' @param x result of enrichment analysis
#'
#' @noRd
has_pairsim <- function(x) {
    if (length(x@termsim) == 0) {
        error_message <- paste(
            "Term similarity matrix not available.",
            "Please use pairwise_termsim function to",
            "deal with the results of enrichment analysis."
        )
        stop(error_message)
    }
}


#' Get graph_from_data_frame() result
#'
#' @importFrom igraph graph.empty
#' @importFrom igraph graph_from_data_frame
#' @param enrichDf A data.frame of enrichment result.
#' @param geneSets A list gene sets with the names of enrichment IDs
#' @param color a string, the column name of y for nodes colours
#' @param cex_line Numeric, scale of line width
#' @param min_edge The minimum similarity threshold for whether
#' two nodes are connected; should be between 0 and 1 (default `0.2`).
#' @param pair_sim Semantic similarity matrix.
#' @param method Method of calculating the similarity between nodes,
#' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
#' "JC" (Jaccard similarity coefficient) methods
#' @return result of graph_from_data_frame()
#' @importFrom igraph V
#' @importFrom igraph 'V<-'
#' @importFrom igraph E
#' @importFrom igraph 'E<-'
#' @importFrom igraph add_vertices
#' @importFrom igraph delete.edges
#' @noRd
build_emap_graph <- function(
    enrichDf,
    geneSets,
    color,
    cex_line,
    min_edge,
    pair_sim
) {
    if (!is.numeric(min_edge) || min_edge < 0 || min_edge > 1) {
        stop('"min_edge" should be a number between 0 and 1.')
    }

    if (is.null(dim(enrichDf)) || nrow(enrichDf) == 1) {
        # when just one node
        g <- igraph::make_empty_graph(n = 1, directed = FALSE)
        V(g)$name <- unname(get_term_labels(enrichDf, enrichDf$ID))
        V(g)$color <- "red"
        return(g)
    } else {
        w <- pair_sim[
            unname(get_term_labels(enrichDf, enrichDf$ID)),
            unname(get_term_labels(enrichDf, enrichDf$ID))
        ]
    }

    wd <- reshape2::melt(w)
    wd <- wd[wd[, 1] != wd[, 2], ]
    # remove NA
    wd <- wd[!is.na(wd[, 3]), ]

    label_map <- get_term_labels(enrichDf, enrichDf$ID)
    keep_edge <- wd[, 3] >= min_edge
    vertex_df <- data.frame(
        name = unname(label_map),
        stringsAsFactors = FALSE
    )
    if (!any(keep_edge)) {
        g <- graph_from_data_frame(
            data.frame(from = character(0), to = character(0)),
            directed = FALSE,
            vertices = vertex_df
        )
    } else {
        edge_df <- wd[keep_edge, , drop = FALSE]
        g <- graph_from_data_frame(edge_df[, -3], directed = FALSE, vertices = vertex_df)
        E(g)$width <- sqrt(edge_df[, 3] * 5) * cex_line
        # Use similarity as the weight(length) of an edge
        E(g)$weight <- edge_df[, 3]
    }

    idx <- match(V(g)$name, label_map)
    cnt <- sapply(geneSets[enrichDf$ID[idx]], length)
    V(g)$size <- cnt
    if (color %in% names(enrichDf)) {
        colVar <- enrichDf[idx, color]
    } else {
        colVar <- color
    }
    V(g)$color <- colVar
    return(g)
}

#' Merge the compareClusterResult file
#'
#' @param yy A data.frame of enrichment result.
#'
#' @return a data.frame
#' @noRd
merge_compareClusterResult <- function(yy) {
    yy_union <- yy[!duplicated(yy$ID), ]
    yy_ids <- lapply(split(yy, yy$ID), function(x) {
        ids <- unique(unlist(strsplit(x$geneID, "/")))
        cnt <- length(ids)
        list(ID = paste0(ids, collapse = "/"), cnt = cnt)
    })

    ids <- vapply(yy_ids, function(x) x$ID, character(1))
    cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))

    yy_union$geneID <- ids[yy_union$ID]
    yy_union$Count <- cnt[yy_union$ID]
    yy_union$Cluster <- NULL
    yy_union
}






#' Get the location of group label
#'
#' @param node_data node information data frame
#' @param label_format A numeric value sets wrap length, alternatively a
#' custom function to format axis labels.
#' @return a data.frame object.
#' @noRd
get_label_location <- function(node_data, label_format) {
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    label_x <- stats::aggregate(x ~ color2, node_data, mean)
    label_y <- stats::aggregate(y ~ color2, node_data, mean)
    data.frame(x = label_x$x, y = label_y$y, label = label_func(label_x$color2))
}


#' Cluster similar nodes together by k-means
#'
#' @param node_data node information data frame.
#' @param enrichDf data.frame of enrichment result.
#' @param nWords Numeric, the number of words in the cluster tags.
#' @param clusterFunction clustering method function, such as `stats::kmeans`, `cluster::clara`,
#' `cluster::fanny`, or `cluster::pam`.
#' @param nCluster Numeric, the number of clusters,
#' the default value is square root of the number of nodes.
#' @noRd
groupNode <- function(
    node_data,
    enrichDf,
    nWords,
    clusterFunction = stats::kmeans,
    nCluster
) {
    wrongMessage <- paste(
        "Wrong clusterFunction parameter or unsupported clustering method;",
        "set to default `clusterFunction = kmeans`"
    )
    if (is.character(clusterFunction)) {
        clusterFunction <- eval(parse(text = clusterFunction))
    }
    if (!"color2" %in% colnames(node_data)) {
        dat <- data.frame(x = node_data$x, y = node_data$y)
        nCluster <- ifelse(
            is.null(nCluster),
            floor(sqrt(nrow(dat))),
            min(nCluster, nrow(dat))
        )
        node_data$color2 <- tryCatch(
            expr = clusterFunction(dat, nCluster)$cluster,
            error = function(e) {
                message(wrongMessage)
                clusterFunction(dat, nCluster)$cluster
            }
        )
        if (is.null(node_data$color2)) {
            message(wrongMessage)
            node_data$color2 <- clusterFunction(dat, nCluster)$cluster
        }
    }
    goid <- enrichDf$ID
    cluster_color <- unique(node_data$color2)
    clusters <- lapply(cluster_color, function(i) {
        goid[which(node_data$color2 == i)]
    })
    cluster_label <- sapply(
        cluster_color,
        get_wordcloud,
        node_data = node_data,
        nWords = nWords
    )
    names(cluster_label) <- cluster_color
    node_data$color2 <- cluster_label[as.character(node_data$color2)]
    return(node_data)   
}

#' Add ellipse to group nodes
#'
#' @param node_data node data frame
#' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
#' The default is FALSE.
#' @param label logical, TRUE to label the ellipse (default)
#' @param ellipse_style style of ellipse, one of "ggforce" and "polygon".
#' @param ellipse_pro numeric indicating confidence value for the ellipses
#' @param alpha the transparency of ellipse fill.
#' @importFrom rlang check_installed
#' @importFrom ggplot2 scale_fill_discrete
#' @noRd
add_ellipse <- function(
    node_data,
    group_legend,
    label = TRUE,
    ellipse_style = "ggforce",
    # ellipse_pro = 0.95,
    alpha = 0.3,
    ...
) {
    show_legend <- c(group_legend, FALSE)
    names(show_legend) <- c("fill", "color")
    ellipse_style <- match.arg(ellipse_style, c("ggforce", "polygon"))

    require_suggested('ggforce', 'for `add_ellipse()`.');

    if (ellipse_style == "ggforce") {
        if (label) {
            p <- ggforce::geom_mark_ellipse(
                data = node_data,
                aes(
                    x = !!sym('x'),
                    y = !!sym('y'),
                    fill = !!sym('color2'),
                    label = !!sym('color2')
                ),
                alpha = alpha,
                color = NA,
                show.legend = show_legend
            )
        } else {
            p <- ggforce::geom_mark_ellipse(
                data = node_data,
                aes(x = !!sym('x'), y = !!sym('y'), fill = !!sym('color2')),
                alpha = alpha,
                color = NA,
                show.legend = show_legend
            )
        }
    }

    # not in used
    if (FALSE && ellipse_style == "polygon") {
        ellipse_pro <- 0.95  # Define default ellipse_pro value
        p <- ggplot2::stat_ellipse(
            data = node_data,
            aes(x = !!sym('x'), y = !!sym('y'), fill = !!sym('color2')),
            geom = "polygon",
            level = ellipse_pro,
            alpha = alpha,
            show.legend = group_legend,
            ...
        )
    }

    if (group_legend) {
        p <- list(p, scale_fill_discrete(name = "groups"))
    }

    return(p)
}


list2df <- ggtangle:::list2df

