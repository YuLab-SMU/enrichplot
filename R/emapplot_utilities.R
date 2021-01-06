##' Get the similarity matrix
##'
##' @param y a data.frame of enrichment result
##' @param geneSets a list, the names of geneSets are term ids,
##' and every object is a vertor of genes
##' @param method method of calculating the similarity between nodes,
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
##' @param y a data.frame of clusterProfiler result
##' @param geneSets a list gene sets with the names of enrichment IDs
##' @param color a string, the column name of y for nodes colours
##' @param cex_line scale of line width
##' @param min_edge minimum percentage of overlap genes to display the edge,
##' should between 0 and 1, default value is 0.2
##' @param pair_sim semantic similarity matrix
##' @param method method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @return result of graph.data.frame()
##' @noRd
build_emap_graph <- function(y, geneSets, color, cex_line, min_edge,
                             pair_sim  = NULL, method = NULL) {

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
##' @param x enrichment result.
##' @param y as.data.frame(x).
##' @param n number of enriched terms to display.
##' @param color variable that used to color enriched terms, e.g. pvalue,
##' p.adjust or qvalue.
##' @param cex_line scale of line width.
##' @param min_edge minimum percentage of overlap genes to display the edge,
##' should between 0 and 1, default value is 0.2.
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
##' @param yy a data.frame of clusterProfiler result
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
##' @param y a data.frame
##' @param g an igraph object
##' @param y_union a data.frame
##' @param cex_category scale of pie plot
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @param layout layout of the map
##' @noRd
build_ggraph <- function(y, g, y_union, cex_category, pie, layout){
    ## when y just have one line
    if(is.null(dim(y)) | nrow(y) == 1) {
        title <- y$Cluster
        p <- ggraph(g, "tree") + geom_node_point(color="red", size=5 * cex_category) +
            geom_node_text(aes_(label=~name)) + theme_void() +
            labs(title=title)
        return(p)
    }

    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        ##return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
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
    ggraph(g, layout=layout)
}

##' Keep selected category in enrichment result
##'
##' @param showCategory a number or a vectory of enriched terms to display
##' @param x enrichment result
##' @param split separate result by 'category' variable
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

##' convert a list of gene IDs to igraph object.
##'
##'
##' @title convert gene IDs to igraph object
##' @param inputList a list of gene IDs
##' @return a igraph object.
##' @importFrom igraph graph.data.frame
##' @author Guangchuang Yu
##' @noRd
list2graph <- function(inputList) {
    x <- list2df(inputList)
    g <- graph.data.frame(x, directed=FALSE)
    return(g)
}

##' convert a list of gene IDs to data.frame object.
##'
##'
##' @title convert gene IDs to data.frame object
##' @param inputList a list of gene IDs
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
##' @param pdata2 data of a ggraph object
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @return a data.frame object.
##' @noRd
get_label_location <- function(pdata2, label_format) {
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    label_x <- stats::aggregate(x ~ color, pdata2, mean)
    label_y <- stats::aggregate(y ~ color, pdata2, mean)
    data.frame(x = label_x$x, y = label_y$y,
        label = label_func(label_x$color))
}

##' Add group label to a ggplot2 object
##'
##' @param repel a logical value, whether to correct the position of the label.
##' @param shadowtext a logical value, whether to use shadow font. 
##' @param p a ggplot2 object.
##' @param label_location a data.frame with the location of group label.
##' @param label_group a numeric value, default size of group label.
##' @param cex_label_group scale of group labels size.
##' @param ... additional parameters.
##' @return a ggplot2 object.
##' @noRd
add_group_label <- function(repel, shadowtext, p, label_location, 
                            label_group, cex_label_group, ...) {
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
##' @param label_location a data.frame with the location of group label.
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





