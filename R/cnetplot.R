##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "enrichResult"),
          function(x, showCategory = 5,
                   foldChange = NULL, layout = "kk", ...) {
              cnetplot.enrichResult(x, showCategory = showCategory,
                                    foldChange = foldChange, layout = layout, ...)
          })

##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "gseaResult"),
          function(x, showCategory = 5,
                   foldChange = NULL, layout = "kk", ...) {
              cnetplot.enrichResult(x, showCategory = showCategory,
                                    foldChange = foldChange, layout = layout, ...)
          })

##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 5,
                   foldChange = NULL, layout = "kk", ...) {
              cnetplot.compareClusterResult(x, showCategory = showCategory,
                                    foldChange = foldChange, layout = layout, ...)
          })


##' @rdname cnetplot
##' @param colorEdge whether coloring edge by enriched terms
##' @param circular whether using circular layout
##' @param node_label select which labels to be displayed.
##' one of 'category', 'gene', 'all' and 'none', default is "all".
##' @param cex_category number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @param cex_gene number indicating the amount by which plotting gene nodes
##' should be scaled relative to the default.
##' @param node_label_size size of node label, this parameter has been
##' changed to cex_label_category and cex_label_gene
##' @param cex_label_category scale of category node label size
##' @param cex_label_gene scale of gene node label size
##' @importFrom ggraph geom_edge_arc
##' @importFrom ggplot2 scale_colour_gradient2
##' @author Guangchuang Yu
cnetplot.enrichResult <- function(x,
                     showCategory = 5,
                     foldChange   = NULL,
                     layout = "kk",
                     colorEdge = FALSE,
                     circular = FALSE,
                     node_label = "all",
                     cex_category = 1,
                     cex_gene = 1,
                     node_label_size = NULL,
                     cex_label_category = 1,
                     cex_label_gene = 1,
                     ...) {

    if (!is.null(node_label_size))
        message("node_label_size parameter has been changed to 'cex_label_category' and 'cex_label_gene'")
    # if (is.null(5 * cex_label_category)) {
        # if (!is.null(node_label_size)) {
            # 5 * cex_label_category <- node_label_size
        # } else {
            # 5 * cex_label_category <- 5
        # }
    # }

    # if (is.null(5 * cex_label_gene)) {
        # if (!is.null(node_label_size)) {
            # 5 * cex_label_gene <- node_label_size
        # } else {
            # 5 * cex_label_gene <- 5
        # }
    # }

    label_category <- 5
    label_gene <- 5
    node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
    if (circular) {
        layout <- "linear"
        geom_edge <- geom_edge_arc
    } else {
        geom_edge <- geom_edge_link
    }

    geneSets <- extract_geneSets(x, showCategory)

    g <- list2graph(geneSets)

    foldChange <- fc_readable(x, foldChange)

    size <- sapply(geneSets, length)
    V(g)$size <- min(size)/2
    n <- length(geneSets)
    V(g)$size[1:n] <- size
    node_scales <- c(rep(cex_category, n), rep(cex_gene, (length(V(g)) - n)))
    if (colorEdge) {
        E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
        edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)
    } else {
        edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
    }

    if (!is.null(foldChange)) {
        fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
        V(g)$color <- NA
        V(g)$color[(n+1):length(V(g))] <- fc
        show_legend <- c(TRUE, FALSE)
        names(show_legend) <- c("color", "size")
        p <- ggraph(g, layout=layout, circular = circular)
        p <- p + edge_layer +
            # geom_node_point(aes_(color=~as.numeric(as.character(color)),
            geom_node_point(aes_(color=~I("#E5C494"), size=~size),
                data = p$data[1:n, ]) +
            scale_size(range=c(3, 8) * cex_category) +
            ggnewscale::new_scale("size") +
            ggnewscale::new_scale_color() +
            geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size),
                data = p$data[-(1:n), ], show.legend = show_legend) +
            scale_size(range=c(3, 3) * cex_gene) +
            scale_colour_gradient2(name = "fold change", low = "blue",
                                   mid = "white", high = "red")
    } else {
        V(g)$color <- "#B3B3B3"
        V(g)$color[1:n] <- "#E5C494"
        p <- ggraph(g, layout=layout, circular=circular)
        p <- p + edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size), data = p$data[1:n, ]) +
            scale_size(range=c(3, 8) * cex_category) +
            ggnewscale::new_scale("size") +
            geom_node_point(aes_(color=~I(color), size=~size),
                data = p$data[-(1:n), ], show.legend = FALSE) +
            scale_size(range=c(3, 3) * cex_gene)
    }

    p <- p + theme_void()

# gan jue zhe li duo chi yi ju, zhi qian yi jing ba bu xu yao de  she zhi cheng le ""
    if (node_label == "category") {
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,],
                size = label_category * cex_label_category, bg.color = "white")
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,],
                size = label_category * cex_label_category)
        }
    } else if (node_label == "gene") {
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                repel=TRUE, size = label_gene * cex_label_gene, bg.color = "white")
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
            repel=TRUE, size = label_gene * cex_label_gene)
        }
    } else if (node_label == "all") {
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                    repel=TRUE, size = label_gene * cex_label_gene, bg.color = "white") + 
                geom_node_text(aes_(label=~name), repel=TRUE,
                    size = label_category * cex_label_category, bg.color = "white", data = p$data[1:n,]) 
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                    repel=TRUE, size = label_gene * cex_label_gene) + 
                geom_node_text(aes_(label=~name), data = p$data[1:n,],
                    repel=TRUE, size = label_category * cex_label_category)
        }

    }

    return(p)
}

##' @param colorEdge whether coloring edge by enriched terms
##' @param circular whether using circular layout
##' @param node_label select which labels to be displayed.
##'                   one of 'category', 'gene', 'all' and 'none', default is "all".
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @param pie_scale scale of pie chart, this parameter has been changed to "node_scale"
##' @param legend_n number of circle in legend
##' @param x_loc,y_loc the location of scatterpie legend
##' @importFrom ggraph geom_edge_arc
##' @noRd
cnetplot.compareClusterResult <- function(x,
                     showCategory = 5,
                     foldChange   = NULL,
                     layout = "kk",
                     colorEdge = FALSE,
                     circular = FALSE,
                     node_label = "all",
                     split=NULL,
                     pie = "equal",
                     pie_scale = NULL,
                     cex_category = 1,
                     cex_gene = 1,
                     legend_n = 5,
                     node_label_size = NULL,
                     x_loc = NULL,
                     y_loc = NULL,
                     cex_label_category = 1,
                     cex_label_gene = 1,
                     ...) {

    if (!is.null(node_label_size))
        message("node_label_size parameter has been changed to 'cex_label_category' and 'cex_label_gene'")
    # if (is.null(5 * cex_label_category)) {
        # if (!is.null(node_label_size)) {
            # cex_label_category <- node_label_size
        # } else {
            # cex_label_category <- 2.5
        # }
    # }

    # if (is.null(5 * cex_label_gene)) {
        # if (!is.null(node_label_size)) {
            # cex_label_gene <- node_label_size
        # } else {
            # cex_label_gene <- 2.5
        # }
    # }

    
    if (!is.null(pie_scale))
        message("pie_scale parameter has been changed to 'cex_category' and 'cex_gene'")
    if (is.null(cex_category)) {
        if (!is.null(pie_scale)) {
            cex_category <- pie_scale
        } else {
            cex_category <- 1
        }
    }

    if (is.null(cex_gene)) {
        if (!is.null(pie_scale)) {
            cex_gene <- pie_scale
        } else {
            cex_gene <- 1
        }
    }
    
    label_category <- 2.5
    label_gene <- 2.5
    range_category_size <- c(3, 8)
    range_gene_size <- c(3, 3)
    y <- fortify(x, showCategory=showCategory,
        includeAll=TRUE, split=split)
    y$Cluster <- sub("\n.*", "", y$Cluster)


    y_union <- get_y_union(y = y, showCategory = showCategory)
    y <- y[y$ID %in% y_union$ID, ]
    node_label <- match.arg(node_label, c("category", "gene", "all", "none"))


    if (circular) {
        layout <- "linear"
        geom_edge <- geom_edge_arc
    } else {
        geom_edge <- geom_edge_link
    }


    #geneSets <- extract_geneSets(x, showCategory)
    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/",
                                  fixed = TRUE), y_union$Description)
    n <- length(geneSets)
    g <- list2graph(geneSets)
    edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
    if(is.null(dim(y)) | nrow(y) == 1) {
        V(g)$size <- 1
        V(g)$size[1] <- 3
        V(g)$color <- "#B3B3B3"
        V(g)$color[1] <- "#E5C494"
        title <- y$Cluster
        p <- ggraph(g, layout=layout, circular=circular)
        p <- p + edge_layer + theme_void() +
            # geom_node_point(aes_(color=~I(color), size=~size)) +
            # labs(title= title) +
            # scale_size(range=c(3, 8) * mean(node_scales)) + theme(legend.position="none")
            geom_node_point(aes_(color=~I(color), size=~size),
                data = p$data[1:n, ]) +
            scale_size(range = range_category_size * cex_category) +
            ggnewscale::new_scale("size") +
            geom_node_point(aes_(color=~I(color), size=~size),
                data = p$data[-(1:n), ], show.legend = FALSE) +
            scale_size(range = range_gene_size * cex_gene) +
            labs(title= title) +
            theme(legend.position="none")
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                    size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
                geom_node_text(aes_(label=~name), data = p$data[1:n,],
                    size = label_category * cex_label_category, bg.color = "white", repel=TRUE)
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                    size = label_gene * cex_label_gene, repel=TRUE) + 
                geom_node_text(aes_(label=~name), data = p$data[1:n,],
                    size = label_category * cex_label_category, repel=TRUE)
        }


        return(p)
    }

    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        p <- ggraph(g) + edge_layer
    } else {
        p <- ggraph(g, layout=layout, circular=circular) + edge_layer
    }


    #pie chart begin
    #obtain the cluster distribution of each GO term and gene
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)

    gene_Cluster_mat <- prepare_pie_gene(y)
    if(ncol(ID_Cluster_mat) > 1) {
        clusters <- match(colnames(ID_Cluster_mat),colnames(gene_Cluster_mat))
        ID_Cluster_mat <- ID_Cluster_mat[,clusters]
        gene_Cluster_mat <- gene_Cluster_mat[,clusters]
    }
    ID_Cluster_mat2 <- rbind(ID_Cluster_mat,gene_Cluster_mat)
    #add the coordinates
    aa <- p$data
    ii <- match(rownames(ID_Cluster_mat2), aa$name)

    ID_Cluster_mat2$x <- aa$x[ii]
    ID_Cluster_mat2$y <- aa$y[ii]
    #add the radius of the pie chart, the radius of go terms mean the number of genes
    ii <- match(rownames(ID_Cluster_mat2)[1:n], y_union$Description)
    node_scales <- c(rep(cex_category, n), rep(cex_gene, (length(V(g)) - n)))
    # sum_yunion <- sum(y_union[ii,9])
    sum_yunion <- sum(y_union[ii, "Count"])
    sizee <- sqrt(y_union[ii, "Count"] / sum_yunion)
    ID_Cluster_mat2$radius <- min(sizee)/2  * sqrt(cex_gene)
    ID_Cluster_mat2$radius[1:n] <- sizee * sqrt(cex_category)
    if(is.null(x_loc)) x_loc <- min(ID_Cluster_mat2$x)
    if(is.null(y_loc)) y_loc <- min(ID_Cluster_mat2$y)
    #node_label
    if (node_label == "category") {
        p$data$name[(n+1):nrow(p$data)] <- ""
    } else if (node_label == "gene") {
        p$data$name[1:n] <- ""
    } else if (node_label == "none") {
        p$data$name <- ""
    }
    if(ncol(ID_Cluster_mat2) > 4) {
    ## should not have foldChange
        if (!is.null(foldChange)) {
            log_fc <- abs(foldChange)
            genes <- rownames(ID_Cluster_mat2)[(n+1):nrow(ID_Cluster_mat2)]
            gene_fc <- rep(1,length(genes))
            gid <- names(log_fc)
            #Turn the id of  gid into gene symbols
            ii <- gid %in% names(x@gene2Symbol)
            gid[ii] <- x@gene2Symbol[gid[ii]]
            ii <- match(genes,gid)
            gene_fc <- log_fc[ii]
            gene_fc[is.na(gene_fc)] <- 1
            gene_fc2 <- c(rep(1,n),gene_fc)
            #Assign value to the size of the genes
            # ID_Cluster_mat2$radius <- min(sizee)/2*gene_fc2
            # ID_Cluster_mat2$radius[1:n] <- sizee

            ID_Cluster_mat2$radius <- min(sizee)/2*gene_fc2 * sqrt(cex_gene)
            ID_Cluster_mat2$radius[1:n] <- sizee * sqrt(cex_category)
            # p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius),
                    # data=ID_Cluster_mat2,
                    # cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)],
                    # color=NA) +
            p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius),
                    data=ID_Cluster_mat2[1:n, ],
                    cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)], color=NA) +
                geom_scatterpie_legend(ID_Cluster_mat2$radius[1:n],
                    x=x_loc, y=y_loc + 3, n = legend_n, labeller=function(x) round(x^2 * sum_yunion / cex_category))  +
                geom_scatterpie(aes_(x=~x,y=~y,r=~radius),
                    data=ID_Cluster_mat2[-(1:n), ],
                    cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)],
                    color=NA, show.legend = FALSE) +
                coord_equal()+
                geom_scatterpie_legend(ID_Cluster_mat2$radius[(n+1):nrow(ID_Cluster_mat2)],
                    x=x_loc, y=y_loc, n = legend_n,
                    labeller=function(x) round(x*2/(min(sizee))/sqrt(cex_gene),3)) +
                ggplot2::annotate("text", x = x_loc + 3, y = y_loc, label = "log2FC")  +
                ggplot2::annotate("text", x = x_loc + 3, y = y_loc + 3, label = "gene number")

            if (utils::packageVersion("ggrepel") >= "0.9.0") {
                p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                        size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
                    geom_node_text(aes_(label=~name), data = p$data[1:n,],
                        size = label_category * cex_label_category, bg.color = "white", repel=TRUE)
            } else {
                p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                        size = label_gene * cex_label_gene, repel=TRUE) + 
                    geom_node_text(aes_(label=~name), data = p$data[1:n,],
                        size = label_category * cex_label_category, repel=TRUE)
            }


            p <- p + theme_void() + labs(fill = "Cluster")
            return(p)
        }
        ## should not have foldChange
        # p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat2,
                # cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)],
                # color=NA) +
           # coord_equal()
        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius),
                data=ID_Cluster_mat2[1:n, ],
                cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)], color=NA) +
            geom_scatterpie(aes_(x=~x,y=~y,r=~radius),
                data=ID_Cluster_mat2[-(1:n), ],
                cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)],
                color=NA, show.legend = FALSE) +
            coord_equal() +
            geom_scatterpie_legend(ID_Cluster_mat2$radius[1:n],
                    x=x_loc, y=y_loc, n = legend_n, labeller=function(x) round(x^2 * sum_yunion / cex_category)) +
            ggplot2::annotate("text", x = x_loc + 3, y = y_loc, label = "gene number")

        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                    size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
                geom_node_text(aes_(label=~name), data = p$data[1:n,],
                    size = label_category * cex_label_category, bg.color = "white", repel=TRUE)
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                    size = label_gene * cex_label_gene, repel=TRUE) + 
                geom_node_text(aes_(label=~name), data = p$data[1:n,],
                    size = label_category * cex_label_category, repel=TRUE)
        }




        p <- p + theme_void() + labs(fill = "Cluster")
        return(p)
    }
    title <- colnames(ID_Cluster_mat2)[1]
    V(g)$size <- ID_Cluster_mat2$radius
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"

    p <- ggraph(g, layout=layout, circular=circular)
    p <- p + edge_layer +
        geom_node_point(aes_(color=~I(color), size=~size),
            data = p$data[1:n, ]) +
        scale_size(range = range_category_size * cex_category) +
        ggnewscale::new_scale("size") +
        geom_node_point(aes_(color=~I(color), size=~size),
            data = p$data[-(1:n), ], show.legend = FALSE) +
        scale_size(range = range_gene_size * cex_gene) +
        labs(title= title)

    if (utils::packageVersion("ggrepel") >= "0.9.0") {
        p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
            geom_node_text(aes_(label=~name), data = p$data[1:n,],
                size = label_category * cex_label_category, bg.color = "white", repel=TRUE)
    } else {
        p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
                size = label_gene * cex_label_gene, repel=TRUE) + 
            geom_node_text(aes_(label=~name), data = p$data[1:n,],
                size = label_category * cex_label_category, repel=TRUE)
    }


    p + theme_void() + theme(legend.position="none")
}




##' convert a list of gene IDs to igraph object.
##'
##'
##' @title convert gene IDs to igraph object
##' @param inputList a list of gene IDs
##' @return a igraph object.
##' @importFrom igraph graph.data.frame
##' @author Guangchuang Yu
list2graph <- function(inputList) {
    x <- list2df(inputList)
    g <- graph.data.frame(x, directed=FALSE)
    return(g)
}


list2df <- function(inputList) {
    # ldf <- lapply(1:length(inputList), function(i) {
    ldf <- lapply(seq_len(length(inputList)), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}


