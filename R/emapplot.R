##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "enrichResult"),
          function(x, showCategory = 30, color = "p.adjust",
                   layout = "nicely", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, color = "p.adjust",
                   layout = "nicely", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 30, color = "p.adjust",
                   layout = "nicely", ...) {

              emapplot.compareClusterResult(x, showCategory = showCategory,
                                            color=color, layout = layout, ...)
          })


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
get_ww <- function(y, geneSets, method, semData = NULL) {
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
    return(w)
}


#' Check whether the similarity matrix exists
#'
#' @param x result of enrichment analysis
#'
has_pairsim <- function(x) {
    if (length(x@termsim) == 0) {
        error_message <- paste("Term similarity matrix not available.",
            "Please use pairwise_termsim function to",
            "deal with the results of enrichment analysis.")
        # error_message <- gsub("[ ]+", " ", error_message)
        # error_message <- gsub("[\r\n]", "", error_message)
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
emap_graph_build <- function(y, geneSets, color, cex_line, min_edge, 
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
        w <- pair_sim
        if (method == "JC") {
            w <- w[as.character(y$Description), as.character(y$Description)]
        } else {
            w <- w[y$ID, y$ID]
        }
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
    
    g <- emap_graph_build(y = y, geneSets = geneSets, color = color,
             cex_line = cex_line, min_edge = min_edge,
             pair_sim = x@termsim, method = x@method)
}



##' @rdname emapplot
##' @importFrom igraph graph.empty
##' @importFrom igraph add_vertices
##' @importFrom igraph graph.data.frame
##' @importFrom igraph delete.edges
##' @importFrom igraph V "V<-"
##' @importFrom igraph E "E<-"
##' @importFrom reshape2 melt
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @param node_scale scale of node, this parameter has been changed to cex_category
##' @param line_scale scale of line width, this parameter has been changed to cex_line
##' @param cex_line scale of line width
##' @param min_edge minimum percentage of overlap genes to display the edge,
##' should between 0 and 1, default value is 0.2
##' @param node_label_size size of node label, this parameter has been
##' changed to cex_label_category
##' @param cex_label_category scale of category node label size
##' @param cex_category number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @author Guangchuang Yu
emapplot.enrichResult <- function(x, showCategory = 30, color="p.adjust",
    layout = "nicely", node_scale = NULL, line_scale = NULL, min_edge=0.2,
    node_label_size = NULL, cex_label_category  = 1, cex_category = NULL,
    cex_line = NULL) {
    
    has_pairsim(x)
    if (!is.null(node_label_size)) 
        message("node_label_size parameter has been changed to 'cex_label_category'")
    # if (is.null(cex_label_category)) {
        # if (!is.null(node_label_size)) {
            # cex_label_category <- node_label_size
        # } else {
            # cex_label_category <- 5
        # }
    # }

    if (!is.null(node_scale)) 
        message("node_scale parameter has been changed to 'cex_category'")
    if (is.null(cex_category)) {
        if (!is.null(node_scale)) {
            cex_category <- node_scale
        } else {
            cex_category <- 1
        }
    }

    if (!is.null(line_scale)) 
        message("line_scale parameter has been changed to 'cex_line'")
    if (is.null(cex_line)) {
        if (!is.null(line_scale)) {
            cex_line <- line_scale
        } else {
            cex_line <- 1
        }
    }
    label_category <- 5
    n <- update_n(x, showCategory)
    # geneSets <- geneInCategory(x) ## use core gene for gsea result


    y <- as.data.frame(x)
    g <- get_igraph(x=x, y=y, n=n, color=color, cex_line=cex_line,
                    min_edge=min_edge)
    if(n == 1) {
        return(ggraph(g,"tree") + geom_node_point(color="red", size=5) +
               geom_node_text(aes_(label=~name)))
    }

    p <- ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }
    p <- p + geom_node_point(aes_(color=~color, size=~size))

    if (utils::packageVersion("ggrepel") >= "0.9.0") {
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
            size = label_category * cex_label_category, bg.color = "white")
    } else {
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
            size = label_category * cex_label_category)
    }
        # geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
    p + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color,
                               guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8) * cex_category)

}



##' Merge the compareClusterResult file
##'
##' @param yy a data.frame of clusterProfiler result
##'
##' @return a data.frame
##' @noRd
merge_compareClusterResult <- function(yy) {
    yy_union<- yy[!duplicated(yy$ID),]
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

##' @rdname emapplot
##' @importFrom igraph E "E<-"
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggplot2 coord_equal
##' @importFrom ggplot2 labs
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom scatterpie geom_scatterpie
##' @importFrom scatterpie geom_scatterpie_legend
##' @importClassesFrom DOSE compareClusterResult
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @param legend_n number of circle in legend
##' @param pie_scale scale of pie chart or point, this parameter has been changed to "node_scale"
##' @param cex_line scale of line width
##' @param min_edge minimum percentage of overlap genes to display the edge, should between 0 and 1, default value is 0.2
##' @importFrom stats setNames
emapplot.compareClusterResult <- function(x, showCategory = 30,
                                          color = "p.adjust", layout = "nicely",
                                          split=NULL, pie = "equal",
                                          legend_n = 5, cex_category = NULL,
                                          pie_scale = NULL, cex_line = 1,
                                          min_edge=0.2, cex_label_category  = 1,
                                          node_label_size = NULL) {
    has_pairsim(x)
    if (!is.null(node_label_size))
        message("node_label_size parameter has been changed to 'cex_label_category'")
    # if (is.null(cex_label_category)) {
        # if (!is.null(node_label_size)) {
            # cex_label_category <- node_label_size
        # } else {
            # cex_label_category <- 3
        # }
    # }

    if (!is.null(pie_scale))
        message("pie_scale parameter has been changed to 'cex_category'")

    if (is.null(cex_category)) {
        if (!is.null(pie_scale)) {
            cex_category <- pie_scale
        } else {
            cex_category <- 1
        }
    }

    label_category <- 3
    ## pretreatment of x, just like dotplot do
    y <- fortify(x, showCategory=showCategory,
                                      includeAll=TRUE, split=split)
    y$Cluster <- sub("\n.*", "", y$Cluster)
    ## geneSets <- geneInCategory(x) ## use core gene for gsea result

    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- get_y_union(y = y, showCategory = showCategory)
    y <- y[y$ID %in% y_union$ID, ]

    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/",
                                  fixed = TRUE), y_union$ID)

    g <- emap_graph_build(y=y_union,geneSets=geneSets,color=color,
                          cex_line=cex_line, min_edge=min_edge,
                          pair_sim = x@termsim, method = x@method)

    p <- get_p(y = y, g = g, y_union = y_union, cex_category = cex_category,
               pie = pie, layout = layout)
    if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(y_union)) | nrow(y_union) == 1)
        return(p)


    p <- ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }

    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- prepare_pie_category(y,pie=pie)


    # plot the edge
    # get the X-coordinate and y-coordinate of pies
    aa <- p$data

    desc <- y_union$Description[match(rownames(ID_Cluster_mat),
                                      y_union$Description)]
    i <- match(desc, aa$name)

    ID_Cluster_mat$x <- aa$x[i]
    ID_Cluster_mat$y <- aa$y[i]

    #Change the radius value to fit the pie plot
    radius <- NULL
    ID_Cluster_mat$radius <- sqrt(aa$size[i] / sum(aa$size) * cex_category)
    #ID_Cluster_mat$radius <- sqrt(aa$size / pi)

    x_loc1 <- min(ID_Cluster_mat$x)
    y_loc1 <- min(ID_Cluster_mat$y)
    ## x_loc2 <- min(ID_Cluster_mat$x)
    ## y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))
    if(ncol(ID_Cluster_mat) > 4) {
        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
            cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
            coord_equal()
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                size = label_category * cex_label_category, bg.color = "white")
        } else {
            p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                size = label_category * cex_label_category)
        }
        p <- p + theme_void() +
            geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
                n = legend_n,
                labeller=function(x) round(sum(aa$size) * x^2 / cex_category)) +
            labs(fill = "Cluster")
        return(p)
    }
    ## annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")
    title <- colnames(ID_Cluster_mat)[1]
    p + geom_node_point(aes_(color=~color, size=~size))
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
            size = label_category * cex_label_category, bg.color = "white")
    } else {
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
            size = label_category * cex_label_category)
    }
    p + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color,
                               guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8) * cex_category)  +labs(title= title)
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
get_p <- function(y, g, y_union, cex_category, pie, layout){
    ## when y just have one line
    if(is.null(dim(y)) | nrow(y) == 1) {
        title <- y$Cluster
        p <- ggraph(g) + geom_node_point(color="red", size=5 * cex_category) +
            geom_node_text(aes_(label=~name)) + theme_void() +
            labs(title=title)
        return(p)
    }

    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        ##return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
        p <- ggraph(g)
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


get_y_union <- function(y, showCategory){
    y_union <- merge_compareClusterResult(y)

    n <- update_n(y_union, showCategory)
    if (is.numeric(n)) {
        y_union <- y_union[1:n,]
    } else {
        y_union <- y_union[match(n, y_union$Description),]
        n <- length(n)
    }
     if (n == 0) {
        stop("no enriched term found...")
    }

   return(y_union)
}







