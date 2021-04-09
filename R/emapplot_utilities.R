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
    n <- nrow(y)
    y_id <- unlist(strsplit(y$ID[1], ":"))[1]
    ## Choose the method to calculate the similarity
    if (method == "JC") {
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$Description
        for (i in seq_len(n-1)) {
            for (j in (i+1):n) {
                w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
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


##' Get graph.data.frame() result
##'
##' @importFrom igraph graph.empty
##' @importFrom igraph graph.data.frame
##' @param y A data.frame of enrichment result.
##' @param geneSets A list gene sets with the names of enrichment IDs
##' @param color a string, the column name of y for nodes colours
##' @param cex_line Numeric, scale of line width
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##' @param pair_sim Semantic similarity matrix.
##' @param method Method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @return result of graph.data.frame()
##' @noRd
build_emap_graph <- function(y, geneSets, color, cex_line, min_edge,
                             pair_sim, method) {

    if (!is.numeric(min_edge) | min_edge < 0 | min_edge > 1) {
    	stop('"min_edge" should be a number between 0 and 1.')
    }

    if (is.null(dim(y)) | nrow(y) == 1) {  # when just one node
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- as.character(y$Description)
        V(g)$color <- "red"
        return(g)
    } else {
        w <- pair_sim[as.character(y$Description), as.character(y$Description)]
    }

    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    # remove NA
    wd <- wd[!is.na(wd[,3]),]
    if (method != "JC") {
        # map id to names
        wd[, 1] <- y[wd[, 1], "Description"]
        wd[, 2] <- y[wd[, 2], "Description"]
    }

    g <- graph.data.frame(wd[, -3], directed=FALSE)
    E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
    # Use similarity as the weight(length) of an edge
    E(g)$weight <- wd[, 3]
    g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    colVar <- y[idx, color]
    V(g)$color <- colVar
    return(g)
}



##' Get an iGraph object
##'
##' @param x Enrichment result.
##' @param y as.data.frame(x).
##' @param n Number of enriched terms to display.
##' @param color variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @param cex_line Scale of line width.
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##'
##' @return an iGraph object
##' @noRd
get_igraph <- function(x, y,  n, color, cex_line, min_edge){
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    if (is.numeric(n)) {
        y <- y[1:n, ]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }

    if (n == 0) {
        stop("no enriched term found...")
    }

    g <- build_emap_graph(y = y, geneSets = geneSets, color = color,
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

##' Get the an ggraph object
##'
##' @importFrom ggplot2 ylim
##' @param x enrichment result.
##' @param y A data.frame of enrichment result.
##' @param y_union A data.frame of enrichment result.
##' @param cex_category Numeric, scale of pie plot.
##' @param pie Proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'.
##' @param layout Layout of the map.
##' @param color a string, the column name of y for nodes colours
##' @param cex_line Numeric, scale of line width
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##' @param pair_sim Semantic similarity matrix.
##' @param method Method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods.
##' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
##' @noRd
build_ggraph <- function(x, y, y_union, cex_category, pie, layout, coords, cex_line,
                        min_edge, pair_sim, method, with_edge){
    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/",
                              fixed = TRUE), y_union$ID) 
                              
    g <- build_emap_graph(y=y_union,geneSets=geneSets,color="p.adjust",
        cex_line=cex_line, min_edge=min_edge, pair_sim = x@termsim, 
        method = x@method)
    ## when y just have one line
    if(is.null(dim(y)) | nrow(y) == 1) {
        title <- y$Cluster
        p <- ggraph(g, "tree") + geom_node_point(color="red", size=5 * cex_category) +
            geom_node_text(aes_(label=~name)) + theme_void() +
            labs(title=title)
        return(p)
    }

    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        p <- ggraph(g, "tree")
        ID_Cluster_mat <- prepare_pie_category(y, pie=pie)

        ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*cex_category)
        colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
            "x", "y", "radius")


        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
                color=NA)+
            xlim(-3,3) + ylim(-3,3) + coord_equal()+
            geom_node_text(aes_(label=~name), repel=TRUE) +
            theme_void()+labs(fill = "Cluster")
        return(p)

    }
    if (!is.null(layout)) {
        p <- ggraph(g, layout=layout)
    } else {
        p <- ggraph(g, layout="nicely")
        if (!is.null(coords)) {
            ggData <- p$data
            ggData$x <- coords$x
            ggData$y <- coords$y
            p$data <- ggData  
        }
    }
    
    
    ## add edges
    if (with_edge & length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }
    return(p)
}

##' Keep selected category in enrichment result
##'
##' @param showCategory A number or a vector of terms. If it is a number, 
##' the first n terms will be displayed. If it is a vector of terms, 
##' the selected terms will be displayed.
##' @param x Enrichment result
##' @param split Separate result by 'category' variable.
##' @noRd
get_selected_category <- function(showCategory, x, split) {
    if (is.numeric(showCategory)) {
        y <- fortify(x, showCategory = showCategory,
                                      includeAll = TRUE, split = split)

    } else {
        y <- as.data.frame(x)
        y <- y[y$Description %in% showCategory, ]
        y <- fortify(y, showCategory=NULL,
                                      includeAll = TRUE, split = split)
    }
    y$Cluster <- sub("\n.*", "", y$Cluster)
    return(y)
}

##' Convert a list of gene IDs to igraph object.
##'
##'
##' @title Convert gene IDs to igraph object
##' @param inputList A list of gene IDs.
##' @return A igraph object.
##' @importFrom igraph graph.data.frame
##' @author Guangchuang Yu
##' @noRd
list2graph <- function(inputList) {
    x <- list2df(inputList)
    g <- graph.data.frame(x, directed=FALSE)
    return(g)
}

##' Convert a list of gene IDs to data.frame object.
##'
##'
##' @title Convert gene IDs to data.frame object
##' @param inputList A list of gene IDs
##' @return a data.frame object.
##' @noRd
list2df <- function(inputList) {
    # ldf <- lapply(1:length(inputList), function(i) {
    ldf <- lapply(seq_len(length(inputList)), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}

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

##' Add group label to a ggplot2 object
##'
##' @param label_style style of group label, one of "shadowtext" and "ggforce".
##' @param repel a logical value, whether to correct the position of the label.
##' @param shadowtext a logical value, whether to use shadow font. 
##' @param p a ggplot2 object.
##' @param label_location a data.frame with the location of group label.
##' @param label_group a numeric value, default size of group label.
##' @param cex_label_group scale of group labels size.
##' @param ... additional parameters.
##' @return a ggplot2 object.
##' @noRd
add_group_label <- function(label_style, repel, shadowtext, p, label_location, 
                            label_group, cex_label_group, ...) {
    if (label_style != "shadowtext") return(p)
    if (!repel) {
        if (shadowtext) {
            p <- p + geom_shadowtext(data = label_location,
                aes_(x =~ x, y =~ y, label =~ label), colour = "black",
                size = label_group * cex_label_group, bg.color = "white", bg.r = 0.1)
        } else {
            p <- p + geom_text(data = label_location,
                aes_(x =~ x, y =~ y, label =~ label), colour = "black",
                size = label_group * cex_label_group)
        }
        
        return(p)
    }

    if (shadowtext) {
        p <- p + ggrepel::geom_text_repel(data = label_location,
            aes_(x =~ x, y =~ y, label =~ label), colour = "black",
            size = label_group * cex_label_group, bg.color = "white", bg.r = 0.1,
            show.legend = FALSE, ...)
    } else {
        p <- p + ggrepel::geom_text_repel(data = label_location,
            aes_(x =~ x, y =~ y, label =~ label), colour = "black",
            size = label_group * cex_label_group, 
            show.legend = FALSE, ...)
    }
    return(p)   
}

##' Add node label to a ggplot2 object
##'
##' @param p a ggplot2 object.
##' @param data it is uesd as the `data` parameter of function `ggraph::geom_node_text`, a data.frame or NULL.
##' @param label_size_node a numeric value to indicate the font size of the node label.
##' @param cex_label_node a numeric value to indicate the scale of node label size.
##' @param shadowtext  a logical value, whether to use shadow font. 
##' @return a ggplot2 object.
##' @noRd
add_node_label <- function(p, data, label_size_node, cex_label_node, shadowtext) {
    if (shadowtext) {
        p <- p + geom_node_text(aes_(label=~name), data = data,
            size = label_size_node * cex_label_node, bg.color = "white", repel=TRUE)
    } else {
        p <- p + geom_node_text(aes_(label=~name), data = data,
            size = label_size_node * cex_label_node, repel=TRUE)
    }
    return(p)
}




##' Cluster similar nodes together by k-means
##' 
##' @param ggData data of a ggraph object.
##' @param y data.frame of enrichment result.
##' @param nWords Numeric, the number of words in the cluster tags.
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' @noRd
groupNode <- function(ggData, y, nWords, nCluster) {
    if (!"color2" %in% colnames(ggData)) {
        dat <- data.frame(x = ggData$x, y = ggData$y)
        if(is.null(nCluster)){
            ggData$color2 <- stats::kmeans(dat, floor(sqrt(nrow(dat))))$cluster
        } else {
            if(nCluster > nrow(dat)) nCluster <- nrow(dat)
                ggData$color2 <- stats::kmeans(dat, nCluster)$cluster
        }
    }
    goid <- y$ID
    cluster_color <- unique(ggData$color2)
    clusters <- lapply(cluster_color, function(i){goid[which(ggData$color2 == i)]})
    cluster_label <- sapply(cluster_color, get_wordcloud, ggData = ggData,
                            nWords=nWords)
    names(cluster_label) <- cluster_color
    ggData$color2 <- cluster_label[as.character(ggData$color2)] 
    return(ggData)
}

##' add ellipse to group the node
##'
##' @param p ggplot2 object
##' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
##' The default is FALSE.
##' @param label_style style of group label, one of "shadowtext" and "ggforce".
##' @param ellipse_style style of ellipse, one of "ggforce" an "polygon".
##' @param ellipse_pro numeric indicating confidence value for the ellipses
##' @param alpha the transparency of ellipse fill.
##' @noRd
add_ellipse <- function(p, group_legend, label_style) {
    show_legend <- c(group_legend, FALSE)
    names(show_legend) <- c("fill", "color") 
    if (label_style == "shadowtext") {
            p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
                         fill =~ color2), show.legend = show_legend)
        } else {
            p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
                         fill =~ color2, label =~ color2), show.legend = show_legend)
        }
    if (group_legend) p <- p + scale_fill_discrete(name = "groups")
    return(p)
}
# add_ellipse <- function(p, group_legend, label_style, 
    # ellipse_style = "ggforce", ellipse_pro = 0.95, alpha = 0.3) {
    # show_legend <- c(group_legend, FALSE)
    # names(show_legend) <- c("fill", "color") 
    # ellipse_style <- match.arg(ellipse_style, c("ggforce", "polygon"))
    # if (ellipse_style == "ggforce") {
        # if (label_style == "shadowtext") {
            # p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
                        # fill =~ color2), alpha = alpha, show.legend = show_legend)
        # } else {
            # p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color2,
                        # fill =~ color2, label =~ color2), alpha = alpha,
                        # show.legend = show_legend)
        # }
        # if (group_legend) p <- p + scale_fill_discrete(name = "groups")
     # } 
    
    # if (ellipse_style == "polygon") {
        # p <- p + ggplot2::stat_ellipse(aes_(x =~ x, y =~ y, fill =~ color2),
                                       # geom = "polygon", level = ellipse_pro,
                                       # alpha = alpha,
                                       # show.legend = group_legend)
    # }

    # return(p)
# }


##' add category nodes 
##'
##' @param p ggplot2 object
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @param color Variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @noRd
add_category_nodes <- function(p, cex_category, color) {
    p + ggnewscale::new_scale_fill() +
        geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ color,
                                    size =~ size)) +
        scale_size_continuous(name = "number of genes",
                              range = c(3, 8) * cex_category) +
        scale_fill_continuous(low = "red", high = "blue", name = color,
                              guide = guide_colorbar(reverse = TRUE)) + 
        theme(legend.title = element_text(size = 10),
                   legend.text  = element_text(size = 10)) +
        theme(panel.background = element_blank()) 
}



##' Get data for pie plot 
##'
##' @param y A data.frame of enrichment result.
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'
##' @param y_union A data.frame of enrichment result.
##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis.
##' @param ggData data of a ggraph object.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @noRd
get_pie_data <- function(y, pie, y_union, cex_pie2axis, ggData, cex_category) {
    ID_Cluster_mat <- prepare_pie_category(y = y, pie=pie) 
    desc <- y_union$Description[match(rownames(ID_Cluster_mat),
                                      y_union$Description)]
    i <- match(desc, ggData$name)
    ID_Cluster_mat$x <- ggData$x[i]
    ID_Cluster_mat$y <- ggData$y[i]
    ID_Cluster_mat$radius <- sqrt(ggData$size[i] / sum(ggData$size) * cex_category * cex_pie2axis)
    return(ID_Cluster_mat)
}

##' Add category node(pie plot) 
##'
##' @param p ggplot2 object.
##' @param ID_Cluster_mat a matrix data for pie plot
##' @param node_label Select which labels to be displayed,
##' one of 'category', 'group', 'all' and 'none'.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis.
##' @param cex_label_category Scale of category node label size.
##' @param shadowtext a logical value, whether to use shadow font.
##' @param legend_n number of circle in legend
##' @param label_size_category Base size of node label.
##' @noRd
add_pie_node <- function(p, ID_Cluster_mat, node_label, 
                         cex_category, cex_pie2axis,
                         cex_label_category,
                         shadowtext, legend_n,
                         label_size_category) {
    color <- NULL
    if(ncol(ID_Cluster_mat) > 4) {
        p <- p + ggnewscale::new_scale_fill() + 
            geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
            cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
            coord_equal() 

        if (node_label == "all" || node_label == "category") {
            p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
                cex_label_node = cex_label_category, shadowtext = shadowtext)
        }
        
        p <- p + theme_void() +
            geom_scatterpie_legend(ID_Cluster_mat$radius, 
                x=min(ID_Cluster_mat$x), y=min(ID_Cluster_mat$y),
                n = legend_n,
                labeller=function(x) round(sum(p$data$size) * x^2 / cex_category/ cex_pie2axis)) +
            labs(fill = "Cluster")
    } else {
        title <- colnames(ID_Cluster_mat)[1]
        p <- p + theme_void() + geom_node_point(aes_(color=~color, size=~size))
        if (node_label == "all" || node_label == "category") {
            p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
                cex_label_node = cex_label_category, shadowtext = shadowtext)
        }
        p <- p + scale_color_continuous(low="red", high="blue", name = color,
                                        guide=guide_colorbar(reverse=TRUE)) +
            scale_size(range=c(3, 8) * cex_category)  +labs(title= title)  
    }  
    return(p)
}
