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
setMethod("cnetplot", signature(x = "list"),
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
##' @param colorEdge Logical, whether coloring edge by enriched terms, the default value is FALSE. 
##' @param circular Logical, whether using circular layout, the default value is FALSE.
##' @param node_label Select which labels to be displayed.
##' one of 'category', 'gene', 'all'(the default) and 'none'.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default, the default value is 1.
##' @param cex_gene Number indicating the amount by which plotting gene nodes
##' should be scaled relative to the default, the default value is 1.
##' @param cex_label_category Scale of category node label size, the 
##' default value is 1.
##' @param cex_label_gene Scale of gene node label size, the default value is 1.
##' @param color_category Color of category node.
##' @param color_gene Color of gene node.
##' @param shadowtext select which node labels to use shadow font,
##' one of 'category', 'gene', 'all' and 'none', default is 'all'.
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
                     cex_label_category = 1,
                     cex_label_gene = 1,
                     color_category = "#E5C494",
                     color_gene = "#B3B3B3",
                     shadowtext = "all",
                     ...) {

    label_size_category <- 5
    label_size_gene <- 5
    node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
    if (circular) {
        layout <- "linear"
        geom_edge <- geom_edge_arc
    } else {
        geom_edge <- geom_edge_link
    }
    if (is.logical(shadowtext)) {
        shadowtext <- ifelse(shadowtext, "all", "none")
    }
    shadowtext_category <- shadowtext_gene <- FALSE
    if (shadowtext == "all") shadowtext_category <- shadowtext_gene <- TRUE
    if (shadowtext == "category") shadowtext_category <- TRUE
    if (shadowtext == "gene") shadowtext_gene <- TRUE
    geneSets <- extract_geneSets(x, showCategory)

    g <- list2graph(geneSets)

    if (!inherits(x,  "list")) {
        foldChange <- fc_readable(x, foldChange)        
    }

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
        # V(g)$color[1:n] <- color_category
        V(g)$color[(n+1):length(V(g))] <- fc
        show_legend <- c(TRUE, FALSE)
        names(show_legend) <- c("color", "size")
        p <- ggraph(g, layout=layout, circular = circular)
        p$data[-(1:n), "size"] <- 3 * cex_gene

        p <- p + edge_layer +
            geom_node_point(aes_(size=~size), color=I(color_category),
                        data = NULL, show.legend = show_legend,
                        alpha = c(rep(1, n), rep(0, nrow(p$data)-n))) +
            ggnewscale::new_scale_color() +
            geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size),
                data = NULL, alpha = c(rep(0, n), rep(1, nrow(p$data)-n))) +
            scale_size(range=c(3, 8) * cex_category) +  
            scale_colour_gradient2(name = "fold change", low = "blue",
                                   mid = "white", high = "red",
                                   guide = guide_colorbar(order = 2))


        # p <- p + edge_layer +
        #     geom_node_point(aes_(color=~I(color_category), size=~size),
        #         data = p$data[1:n, ]) +
        #     scale_size(range=c(3, 8) * cex_category) +
        #     ggnewscale::new_scale_color() +
        #     geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~I(3 * cex_gene)),
        #                     data = p$data[-(1:n), ], show.legend = show_legend) +
        #     scale_colour_gradient2(name = "fold change", low = "blue",
        #                            mid = "white", high = "red",
        #                            guide = guide_colorbar(order = 2))
    } else {
        V(g)$color <- color_gene
        V(g)$color[1:n] <- color_category
        p <- ggraph(g, layout=layout, circular=circular)
        # p <- p + edge_layer +
        #     geom_node_point(aes_(color=~I(color), size=~size), data = p$data[1:n, ]) +
        #     scale_size(range=c(3, 8) * cex_category) +
        #     geom_node_point(aes_(color=~I(color), size=~I(3 * cex_gene)),
        #                     data = p$data[-(1:n), ], show.legend = FALSE) 
        p$data[-(1:n), "size"] <- 3 * cex_gene
        p <- p + edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size))+
            scale_size(range=c(3, 8) * cex_category) 

    }

    p <- p + theme_void()

    if (node_label == "category") {       
        p$data[1:n, "name"] <- NA     
        p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category)
    } else if (node_label == "gene") {
        p$data[-c(1:n), "name"] <- NA
        p <- add_node_label(p = p, data = NULL, label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
    } else if (node_label == "all") {
        p <- add_node_label(p = p, data = NULL,
            label_size_node = c(rep(label_size_category, n), rep(label_size_gene, nrow(p$data)-n)),
            cex_label_node = c(rep(cex_label_category, n), rep(cex_label_gene,, nrow(p$data)-n)), 
            shadowtext = shadowtext_gene)
        # p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
        #     cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
        # p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
        #     cex_label_node = cex_label_category, shadowtext = shadowtext_category)
    }
    if (!is.null(foldChange)) {
        p <- p + guides(size  = guide_legend(order = 1), 
                        color = guide_colorbar(order = 2))
    }
    return(p)
}


##' @param split Separate result by 'category' variable.
##' @param pie Proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'.
##' @param legend_n Number of circle in legend, the default value is 5.
##' @param x_loc,y_loc The location of scatterpie legend.
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
                     cex_category = 1,
                     cex_gene = 1,
                     legend_n = 5,
                     x_loc = NULL,
                     y_loc = NULL,
                     cex_label_category = 1,
                     cex_label_gene = 1,
                     shadowtext = "all",
                     ...) {

    label_size_category <- 2.5
    label_size_gene <- 2.5
    range_category_size <- c(3, 8)
    range_gene_size <- c(3, 3)
    if (is.logical(shadowtext)) {
        shadowtext <- ifelse(shadowtext, "all", "none")
    }
    shadowtext_category <- shadowtext_gene <- FALSE
    if (shadowtext == "all") shadowtext_category <- shadowtext_gene <- TRUE
    if (shadowtext == "category") shadowtext_category <- TRUE
    if (shadowtext == "gene") shadowtext_gene <- TRUE
    ## If showCategory is a number, keep only the first showCategory of each group,
    ## otherwise keep the total showCategory rows
    # y <- get_selected_category(showCategory, x, split)  
    y <- fortify(x, showCategory = showCategory,
                 includeAll = TRUE, split = split)
    y$Cluster <- sub("\n.*", "", y$Cluster)
    
    if ("core_enrichment" %in% colnames(y)) { ## for GSEA result
        y$geneID <- y$core_enrichment
    }
    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- merge_compareClusterResult(y)
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
            geom_node_point(aes_(color=~I(color), size=~size),
                data = p$data[1:n, ]) +
            scale_size(range = range_category_size * cex_category) +
            ggnewscale::new_scale("size") +
            geom_node_point(aes_(color=~I(color), size=~size),
                data = p$data[-(1:n), ], show.legend = FALSE) +
            scale_size(range = range_gene_size * cex_gene) +
            labs(title= title) +
            theme(legend.position="none")    
        p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
        p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category)   
            # geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
            #     size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
            # geom_node_text(aes_(label=~name), data = p$data[1:n,],
            #     size = label_category * cex_label_category, bg.color = "white", repel=TRUE)
        return(p)
    }

    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        p <- ggraph(g, "tree") + edge_layer
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

            # p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
            #         size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
            #     geom_node_text(aes_(label=~name), data = p$data[1:n,],
            #         size = label_category * cex_label_category, bg.color = "white", repel=TRUE) + 
            p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
                cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
            p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
                cex_label_node = cex_label_category, shadowtext = shadowtext_category)
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
        ## add node label
        # p <- p + geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
        #         size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
        #     geom_node_text(aes_(label=~name), data = p$data[1:n,],
        #         size = label_category * cex_label_category, bg.color = "white", repel=TRUE) + 
        p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
        p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category)
        p <- p + theme_void() + labs(fill = "Cluster")
        return(p)
    }
    title <- colnames(ID_Cluster_mat2)[1]
    V(g)$size <- ID_Cluster_mat2$radius
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"

    ggraph(g, layout=layout, circular=circular) + 
        edge_layer +
        geom_node_point(aes_(color=~I(color), size=~size),
            data = p$data[1:n, ]) +
        scale_size(range = range_category_size * cex_category) +
        ggnewscale::new_scale("size") +
        geom_node_point(aes_(color=~I(color), size=~size),
            data = p$data[-(1:n), ], show.legend = FALSE) +
        scale_size(range = range_gene_size * cex_gene) +
        labs(title= title)
        # geom_node_text(aes_(label=~name), data = p$data[-(1:n),],
        #     size = label_gene * cex_label_gene, bg.color = "white", repel=TRUE) + 
        # geom_node_text(aes_(label=~name), data = p$data[1:n,],
        #     size = label_category * cex_label_category, bg.color = "white", repel=TRUE) + 
        p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
        p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category)
        p <- p + theme_void() + theme(legend.position="none")
}

