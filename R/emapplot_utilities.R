##' Get the similarity matrix
##'
##' @param y A data.frame of enrichment result
##' @param geneSets A list, the names of geneSets are term ids,
##' and every object is a vertor of genes.
##' @param method Method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @param semData GOSemSimDATA object
##' @noRd
get_similarity_matrix <- function(y, geneSets, method, semData = NULL) {
    id <- y[, "ID"]
    geneSets <- geneSets[id]
    y_id <- unlist(strsplit(y$ID[1], ":"))[1]
    ## Choose the method to calculate the similarity
    if (method == "JC") {
        w <- .cal_jc_similarity(geneSets, id = id, name = y$Description)
        return(w)
    }

    if (y_id == "GO") {
        if(is.null(semData)) {
            stop("The semData parameter is missing,
                and it can be obtained through godata function in GOSemSim package.")
        }
        w <- GOSemSim::mgoSim(id, id, semData=semData, measure=method,
                              combine=NULL)
    }

    if (y_id == "DOID") w <- DOSE::doSim(id, id, measure=method)
    rownames(y) <- y$ID
    rownames(w) <- colnames(w) <- y[colnames(w), "Description"]
    return(w)
}


##' Check whether the similarity matrix exists
##'
##' @param x result of enrichment analysis
##'
##' @noRd
has_pairsim <- function(x) {
    if (length(x@termsim) == 0) {
        error_message <- paste("Term similarity matrix not available.",
            "Please use pairwise_termsim function to",
            "deal with the results of enrichment analysis.")
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
#' two nodes are connected, should between 0 and 1, default value is 0.2.
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
build_emap_graph <- function(enrichDf, geneSets, color, cex_line, min_edge,
                             pair_sim, method) {

    if (!is.numeric(min_edge) | min_edge < 0 | min_edge > 1) {
    	stop('"min_edge" should be a number between 0 and 1.')
    }

    if (is.null(dim(enrichDf)) | nrow(enrichDf) == 1) {  # when just one node
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- as.character(enrichDf$Description)
        V(g)$color <- "red"
        return(g)
    } else {
        w <- pair_sim[as.character(enrichDf$Description), 
            as.character(enrichDf$Description)]
    }

    wd <- reshape2::melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    # remove NA
    wd <- wd[!is.na(wd[,3]),]
    if (method != "JC") {
        # map id to names
        wd[, 1] <- enrichDf[wd[, 1], "Description"]
        wd[, 2] <- enrichDf[wd[, 2], "Description"]
    }

    g <- graph_from_data_frame(wd[, -3], directed=FALSE)
    E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
    # Use similarity as the weight(length) of an edge
    E(g)$weight <- wd[, 3]
    g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == enrichDf$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    colVar <- enrichDf[idx, color]
    V(g)$color <- colVar
    return(g)
}



##' Get an iGraph object
##'
##' @param x Enrichment result.
##' @param nCategory Number of enriched terms to display.
##' @param color variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @param cex_line Scale of line width.
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##'
##' @return an iGraph object
##' @noRd
get_igraph <- function(x, nCategory, color, cex_line, min_edge){
    y <- as.data.frame(x)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    if (is.numeric(nCategory)) {
        y <- y[1:nCategory, ]
    } else {
        y <- y[match(nCategory, y$Description),]
        nCategory <- length(nCategory)
    }

    if (nCategory == 0) {
        stop("no enriched term found...")
    }

    g <- build_emap_graph(enrichDf = y, geneSets = geneSets, color = color,
             cex_line = cex_line, min_edge = min_edge,
             pair_sim = x@termsim, method = x@method)
}


##' Merge the compareClusterResult file
##'
##' @param yy A data.frame of enrichment result.
##'
##' @return a data.frame
##' @noRd
merge_compareClusterResult <- function(yy) {
    yy_union <- yy[!duplicated(yy$ID),]
    yy_ids <- lapply(split(yy, yy$ID), function(x) {
        ids <- unique(unlist(strsplit(x$geneID, "/")))
        cnt <- length(ids)
        list(ID=paste0(ids, collapse="/"), cnt=cnt)
    })

    ids <- vapply(yy_ids, function(x) x$ID, character(1))
    cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))

    yy_union$geneID <- ids[yy_union$ID]
    yy_union$Count <- cnt[yy_union$ID]
    yy_union$Cluster <- NULL
    yy_union
}

##' add alpha attribute to ggraph edges
##'
##' @param g ggraph object.
##' @param hilight_category category nodes to be highlight.
##' @param alpha_hilight alpha of highlighted nodes.
##' @param alpha_nohilight alpha of unhighlighted nodes.
##' @noRd
edge_add_alpha <- function(g, hilight_category, alpha_nohilight, alpha_hilight) {
    if (!is.null(hilight_category) && length(hilight_category) > 0) {
        edges <- attr(E(g), "vnames")
        E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
        hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
        E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
        # E(g)$alpha[hilight_edge] <- alpha_hilight
        } else {
            E(g)$alpha <- rep(min(0.8, alpha_hilight), length(E(g)))
    }
    return(g)
}

##' add alpha attribute to ggraph nodes
##'
##' @param p ggraph object.
##' @param hilight_category category nodes to be highlight.
##' @param hilight_gene gene nodes to be highlight.
##' @param alpha_hilight alpha of highlighted nodes.
##' @param alpha_nohilight alpha of unhighlighted nodes.
##' @noRd
node_add_alpha <- function(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight) {
    alpha_node <- rep(1, nrow(p$data))
    if (!is.null(hilight_category)) {
        alpha_node <- rep(alpha_nohilight, nrow(p$data)) 
        hilight_node <- c(hilight_category, hilight_gene)
        alpha_node[match(hilight_node, p$data$name)] <- alpha_hilight
    }
    p$data$alpha <- alpha_node
    return(p)
}


# ' Get the an ggraph object
# '
# ' @importFrom ggplot2 ylim
# ' @param x enrichment result.
# ' @param enrichDf A data.frame of enrichment result.
# ' @param mergedEnrichDf A data.frame of merged enrichment result.
# ' @param cex_category Numeric, scale of pie plot.
# ' @param pie Proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'.
# ' @param layout Layout of the map.
# ' @param color a string, the column name of y for nodes colours
# ' @param cex_line Numeric, scale of line width
# ' @param min_edge The minimum similarity threshold for whether 
# ' two nodes are connected, should between 0 and 1, default value is 0.2.
# ' @param pair_sim Semantic similarity matrix.
# ' @param method Method of calculating the similarity between nodes,
# ' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
# ' "JC" (Jaccard similarity coefficient) methods.
# ' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
# ' @param hilight_category category nodes to be highlight.
# ' @param alpha_hilight alpha of highlighted nodes.
# ' @param alpha_nohilight alpha of unhighlighted nodes.
# ' @noRd
# build_ggraph <- function(x, enrichDf, mergedEnrichDf, cex_category, pie, 
#                          layout, coords, cex_line, min_edge, pair_sim,
#                          method, with_edge, hilight_category, alpha_hilight,
#                          alpha_nohilight){
    
#     segment.size <- get_ggrepel_segsize()
#     geneSets <- setNames(strsplit(as.character(mergedEnrichDf$geneID), "/",
#                               fixed = TRUE), mergedEnrichDf$ID) 
                              
#     g <- build_emap_graph(enrichDf=mergedEnrichDf,geneSets=geneSets,color="p.adjust",
#         cex_line=cex_line, min_edge=min_edge, pair_sim = x@termsim, 
#         method = x@method)
#     g <- edge_add_alpha(g, hilight_category, alpha_nohilight, alpha_hilight)
#     # if (!is.null(hilight_category) && length(hilight_category) > 0) {     
#     #     edges <- attr(E(g), "vnames")
#     #     E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
#     #     hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
#     #     E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
#     #     # E(g)$alpha[hilight_edge] <- alpha_hilight
#     # } else {
#     #     E(g)$alpha <- rep(min(0.8, alpha_hilight), length(E(g)))
#     # }
#     ## when enrichDf just have one line
#     if(is.null(dim(enrichDf)) | nrow(enrichDf) == 1) {
#         title <- enrichDf$Cluster
#         p <- ggraph(g, "tree") + geom_node_point(color="red", size=5 * cex_category) +
#             geom_node_text(aes_(label=~name)) + theme_void() +
#             labs(title=title)
#         return(p)
#     }

#     if(is.null(dim(mergedEnrichDf)) | nrow(mergedEnrichDf) == 1) {
#         p <- ggraph(g, "tree")
#         ID_Cluster_mat <- prepare_pie_category(enrichDf = enrichDf, pie=pie)

#         ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*cex_category)
#         colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
#             "x", "y", "radius")


#         p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
#                 cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
#                 color=NA)+
#             xlim(-3,3) + ylim(-3,3) + coord_equal()+
#             geom_node_text(aes_(label=~name), repel=TRUE, segment.size = segment.size) +
#             theme_void()+labs(fill = "Cluster")
#         return(p)

#     }

#     p <- adj_layout(g = g, layout = layout, coords = coords)    
#     p$data$alpha <- rep(1, nrow(p$data))
#     if (!is.null(hilight_category) && length(hilight_category) > 0) {
#         p$data$alpha <- rep(alpha_nohilight, nrow(p$data))
#         ii <- match(hilight_category, p$data$name)
#         p$data$alpha[ii] <- alpha_hilight      
#     }
#     ## add edges
#     if (with_edge & length(E(g)$width) > 0) {
#         p <- p + geom_edge_link(aes_(width=~I(width), alpha=~I(alpha)),
#                                 colour='darkgrey', show.legend = FALSE)
#     }
#     return(p)
# }


# ' Adjust the layout
# '
# ' @param g ggraph object.
# ' @param layout Layout of the map.
# ' @param coords a data.frame with two columns: 'x' for X-axis coordinate and 'y' for Y-axis coordinate.
# ' @noRd
# adj_layout <- function(g, layout, coords) {
#     if (!is.null(layout)) {
#         p <- ggraph(g, layout=layout)
#     } else {
#         p <- ggraph(g, layout="nicely")
#         if (!is.null(coords)) {
#             # ggData <- p$data
#             # rownames(ggData) <- ggData$name
#             # ids <- intersect(ggData$name, rownames(coords))
#             # if (length(ids) > 0) {
#             #     coords <- coords[ids, ]
#             #     ggData[ids, "x"] <- coords$x
#             #     ggData[ids, "y"] <- coords$y
#             #     rownames(ggData) <- rownames(p$data)
#             #     p$data <- ggData 
#             ids <- match(p$data$name, rownames(coords))
#             ids <- ids[!is.na(ids)]
#             if (length(ids) > 0) {
#                 p$data[ids, "x"] <- coords$x
#                 p$data[ids, "y"] <- coords$y
#             } else {
#                 wm <-  paste("Invalid coords parameter, the rownames of the coords", 
#                     "must be the term names found in the emapplot diagram")
#                 warning(wm)
#             }         
#         }
#     }
#     return(p)
# }




##' Get the location of group label
##'
##' @param ggData data of a ggraph object
##' @param label_format A numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @return a data.frame object.
##' @noRd
get_label_location <- function(ggData, label_format) {
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    label_x <- stats::aggregate(x ~ color2, ggData, mean)
    label_y <- stats::aggregate(y ~ color2, ggData, mean)
    data.frame(x = label_x$x, y = label_y$y,
        label = label_func(label_x$color2))
}

# ##' Add group label to a ggplot2 object
# ##'
# ##' @param label_style style of group label, one of "shadowtext" and "ggforce".
# ##' @param repel a logical value, whether to correct the position of the label.
# ##' @param shadowtext a logical value, whether to use shadow font. 
# ##' @param p a ggplot2 object.
# ##' @param label_location a data.frame with the location of group label.
# ##' @param label_group a numeric value, default size of group label.
# ##' @param cex_label_group scale of group labels size.
# ##' @param ... additional parameters.
# ##' @importFrom rlang check_installed
# ##' @return a ggplot2 object.
# ##' @noRd
# add_group_label <- function(label_style, repel, shadowtext, p, label_location, 
#                             label_group, cex_label_group, ...) {
#     if (label_style != "shadowtext") return(p)
#     segment.size <- get_ggrepel_segsize()
#     if (!repel) {
#         if (shadowtext) {
#             p <- p + geom_shadowtext(data = label_location,
#                 aes_(x =~ x, y =~ y, label =~ label), colour = "black",
#                 size = label_group * cex_label_group, bg.color = "white", bg.r = 0.1)
#         } else {
#             p <- p + geom_text(data = label_location,
#                 aes_(x =~ x, y =~ y, label =~ label), colour = "black",
#                 size = label_group * cex_label_group)
#         }
        
#         return(p)
#     }

#     check_installed('ggrepel', 'for `add_group_label()` with `repel = TRUE`.')  
#     if (shadowtext) {
#         p <- p + ggrepel::geom_text_repel(data = label_location,
#             aes_(x =~ x, y =~ y, label =~ label), colour = "black",
#             size = label_group * cex_label_group, bg.color = "white", bg.r = 0.1,
#             show.legend = FALSE, segment.size = segment.size, ...)
#     } else {
#         p <- p + ggrepel::geom_text_repel(data = label_location,
#             aes_(x =~ x, y =~ y, label =~ label), colour = "black",
#             size = label_group * cex_label_group, 
#             show.legend = FALSE, segment.size = segment.size, ...)
#     }
#     return(p)   
# }

# ##' Add node label to a ggplot2 object
# ##'
# ##' @param p a ggplot2 object.
# ##' @param data it is uesd as the `data` parameter of function `ggraph::geom_node_text`, a data.frame or NULL.
# ##' @param label_size_node a numeric value to indicate the font size of the node label.
# ##' @param cex_label_node a numeric value to indicate the scale of node label size.
# ##' @param shadowtext  a logical value, whether to use shadow font. 
# ##' @return a ggplot2 object.
# ##' @noRd
# add_node_label <- function(p, data, label_size_node, cex_label_node, shadowtext) {
#     # If use 'aes_(alpha =~I(alpha))' will put an error for AsIs object.
#     # But I(alpha) is necessory, so use 'alpha = I(data$alpha)'.
#     segment.size <- get_ggrepel_segsize()
#     if (is.null(data)) {
#         data <- p$data
#     }
#     if (shadowtext) {
#         p <- p + geom_node_text(aes_(label=~name), data = data,
#             alpha = I(data$alpha),
#             size = label_size_node * cex_label_node, bg.color = "white", 
#             repel=TRUE, segment.size = segment.size)
#     } else {
#         p <- p + geom_node_text(aes_(label=~name), data = data,
#             alpha = I(data$alpha),
#             size = label_size_node * cex_label_node, repel=TRUE,
#             segment.size = segment.size)
#     }
#     return(p)
# }

##' Cluster similar nodes together by k-means
##' 
##' @param p a ggraph object.
##' @param enrichDf data.frame of enrichment result.
##' @param nWords Numeric, the number of words in the cluster tags.
##' @param clusterFunction function of Clustering method, such as stats::kmeans, cluster::clara,
##' cluster::fanny or cluster::pam.
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' @noRd
groupNode <- function(p, enrichDf, nWords, clusterFunction =  stats::kmeans, nCluster) {
    ggData <- p$data
    wrongMessage <- paste("Wrong clusterFunction parameter or unsupported clustering method;",
         "set to default `clusterFunction = kmeans`")
    if (is.character(clusterFunction)) {
        clusterFunction <- eval(parse(text=clusterFunction))
    }   
    if (!"color2" %in% colnames(ggData)) {
        dat <- data.frame(x = ggData$x, y = ggData$y)
        nCluster <- ifelse(is.null(nCluster), floor(sqrt(nrow(dat))), 
                           min(nCluster, nrow(dat)))
        ggData$color2 <- tryCatch(expr  = clusterFunction(dat, nCluster)$cluster, 
                                  error = function(e) {
                                      message(wrongMessage)
                                      stats::kmeans(dat, nCluster)$cluster
                                  })
        if (is.null(ggData$color2)) {
            message(wrongMessage)
            ggData$color2 <- stats::kmeans(dat, nCluster)$cluster
        }
    }
    goid <- enrichDf$ID
    cluster_color <- unique(ggData$color2)
    clusters <- lapply(cluster_color, function(i){goid[which(ggData$color2 == i)]})
    cluster_label <- sapply(cluster_color, get_wordcloud, ggData = ggData,
                            nWords=nWords)
    names(cluster_label) <- cluster_color
    ggData$color2 <- cluster_label[as.character(ggData$color2)] 
    return(ggData)
}

#' add ellipse to group the node
#'
#' @param p ggplot2 object
#' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
#' The default is FALSE.
#' @param label_style style of group label, one of "shadowtext" and "ggforce".
#' @param ellipse_style style of ellipse, one of "ggforce" an "polygon".
#' @param ellipse_pro numeric indicating confidence value for the ellipses
#' @param alpha the transparency of ellipse fill.
#' @importFrom rlang check_installed
#' @importFrom ggplot2 scale_fill_discrete
#' @noRd
add_ellipse <- function(p, group_legend, label_style, 
    ellipse_style = "ggforce", ellipse_pro = 0.95, alpha = 0.3, ...) {
    show_legend <- c(group_legend, FALSE)
    names(show_legend) <- c("fill", "color") 
    ellipse_style <- match.arg(ellipse_style, c("ggforce", "polygon"))
    
    check_installed('ggforce', 'for `add_ellipse()`.')
    
    if (ellipse_style == "ggforce") {
        if (label_style == "shadowtext") {
            p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
                        fill =~ color2), alpha = alpha, show.legend = show_legend)
        } else {
            p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
                        fill =~ color2, label =~ color2), alpha = alpha,
                        show.legend = show_legend)
        }
        if (group_legend) p <- p + scale_fill_discrete(name = "groups")
     } 
    
    if (ellipse_style == "polygon") {
        p <- p + ggplot2::stat_ellipse(aes_(x =~ x, y =~ y, fill =~ color2),
                                       geom = "polygon", level = ellipse_pro,
                                       alpha = alpha,
                                       show.legend = group_legend, ...)
    }

    return(p)
}
# add_ellipse <- function(p, group_legend, label_style) {
#     show_legend <- c(group_legend, FALSE)
#     names(show_legend) <- c("fill", "color") 
#     if (label_style == "shadowtext") {
#             p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
#                          fill =~ color2), show.legend = show_legend)
#         } else {
#             p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
#                          fill =~ color2, label =~ color2), show.legend = show_legend)
#         }
#     if (group_legend) p <- p + scale_fill_discrete(name = "groups")
#     return(p)
# }


# ##' add category nodes 
# ##'
# ##' @param p ggplot2 object
# ##' @param cex_category Number indicating the amount by which plotting category
# ##' nodes should be scaled relative to the default.
# ##' @param color Variable that used to color enriched terms, e.g. 'pvalue',
# ##' 'p.adjust' or 'qvalue'.
# ##' @noRd
# add_category_nodes <- function(p, cex_category, color) {
#     p + ggnewscale::new_scale_fill() +
#         geom_point(shape = 21, aes_(x=~x, y=~y, fill=~color,
#                                     size=~size, alpha=~I(alpha))) +
#         scale_size_continuous(name = "number of genes",
#                               range = c(3, 8) * cex_category) +
#         # scale_fill_continuous(name = color) + 
#         set_enrichplot_color(type = "fill", name = color) +
#         theme(legend.title = element_text(size = 10),
#                    legend.text  = element_text(size = 10)) +
#         theme(panel.background = element_blank()) 
# }



##' Get data for pie plot 
##'
##' @param enrichDf A data.frame of enrichment result.
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'
##' @param mergedEnrichDf A data.frame of merged enrichment result.
##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis.
##' @param p a ggraph object.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @noRd
get_pie_data <- function(enrichDf, pie, mergedEnrichDf, cex_pie2axis, p, cex_category) {
    ggData <- p$data
    ID_Cluster_mat <- prepare_pie_category(enrichDf = enrichDf, pie=pie) 
    desc <- mergedEnrichDf$Description[match(rownames(ID_Cluster_mat),
                                      mergedEnrichDf$Description)]
    i <- match(desc, ggData$name)
    ID_Cluster_mat$x <- ggData$x[i]
    ID_Cluster_mat$y <- ggData$y[i]
    ID_Cluster_mat$radius <- sqrt(ggData$size[i] / sum(ggData$size) * cex_category * cex_pie2axis)
    return(ID_Cluster_mat)
}

# ##' Add category node(pie plot) 
# ##'
# ##' @param p ggplot2 object.
# ##' @param ID_Cluster_mat a matrix data for pie plot
# ##' @param node_label Select which labels to be displayed,
# ##' one of 'category', 'group', 'all' and 'none'.
# ##' @param cex_category Number indicating the amount by which plotting category
# ##' nodes should be scaled relative to the default.
# ##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis.
# ##' @param cex_label_category Scale of category node label size.
# ##' @param shadowtext a logical value, whether to use shadow font.
# ##' @param legend_n number of circle in legend
# ##' @param label_size_category Base size of node label.
# ##' @noRd
# add_pie_node <- function(p, ID_Cluster_mat, node_label, 
#                          cex_category, cex_pie2axis,
#                          cex_label_category,
#                          shadowtext, legend_n,
#                          label_size_category) {
#     color <- NULL
#     if(ncol(ID_Cluster_mat) > 4) {
#         ID_Cluster_mat$alpha <- p$data$alpha
#         p <- p + ggnewscale::new_scale_fill() + 
#             geom_scatterpie(aes_(x=~x,y=~y,r=~radius,alpha=~I(alpha)), data=ID_Cluster_mat,
#             cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-4)],color=NA) +
#             coord_equal() 

#         if (node_label == "all" || node_label == "category") {
#             p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
#                 cex_label_node = cex_label_category, shadowtext = shadowtext)
#         }
        
#         p <- p + theme_void() +
#             geom_scatterpie_legend(ID_Cluster_mat$radius, 
#                 x=min(ID_Cluster_mat$x), y=min(ID_Cluster_mat$y),
#                 n = legend_n,
#                 labeller=function(x) round(sum(p$data$size) * x^2 / cex_category/ cex_pie2axis)) +
#             labs(fill = "Cluster")
#     } else {
#         title <- colnames(ID_Cluster_mat)[1]
#         p <- p + theme_void() + geom_node_point(aes_(color=~color, size=~size))
#         if (node_label == "all" || node_label == "category") {
#             p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
#                 cex_label_node = cex_label_category, shadowtext = shadowtext)
#         }
#         p <- p + # scale_color_continuous(name = color) +
#             set_enrichplot_color(name = color) + 
#             scale_size(range=c(3, 8) * cex_category)  +labs(title= title)  
#     }  
#     return(p)
# }

list2df <- ggtangle:::list2df
