##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "enrichResult"),
          function(x,  ...) {
              cnetplot.enrichResult(x,  ...)
          })

##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "list"),
          function(x, ...) {
              cnetplot.enrichResult(x, ...)
          })

##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "gseaResult"),
          function(x, ...) {
              cnetplot.enrichResult(x, ...)
          })

##' @rdname cnetplot
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x = "compareClusterResult"),
          function(x, ...) {
              cnetplot.compareClusterResult(x, ...)
          })


##' @rdname cnetplot
##' @param colorEdge Logical, whether coloring edge by enriched terms, the default value is FALSE. 
##' Will be removed in the next version.
##' @param circular Logical, whether using circular layout, the default value is FALSE.
##' Will be removed in the next version.
##' @param node_label Select which labels to be displayed.
##' one of 'category', 'gene', 'all'(the default) and 'none'.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default, the default value is 1.
##' Will be removed in the next version.
##' @param cex_gene Number indicating the amount by which plotting gene nodes
##' should be scaled relative to the default, the default value is 1.
##' Will be removed in the next version.
##' @param cex_label_category Scale of category node label size, the 
##' default value is 1.
##' Will be removed in the next version.
##' @param cex_label_gene Scale of gene node label size, the default value is 1.
##' Will be removed in the next version.
##' @param color_category Color of category node.
##' Will be removed in the next version.
##' @param color_gene Color of gene node.
##' Will be removed in the next version.
##' @param shadowtext select which node labels to use shadow font,
##' one of 'category', 'gene', 'all' and 'none', default is 'all'.

##' @param color.params list, the parameters to control the attributes of highlighted nodes and edges.
##' see the color.params in the following.
##' color.params control the attributes of highlight, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{foldChange} Fold Change of nodes for enrichResult, or size of nodes for compareClusterResult, 
##'         the default value is NULL. 
##'         \item \code{edge} Logical, whether coloring edge by enriched terms, the default value is FALSE. 
##'         \item \code{category} Color of category node.
##'         \item \code{gene} Color of gene node.
##'     }
##' @param cex.params list, the parameters to control the size of nodes and lables.
##' see the cex.params in the following.
##' cex.params control the attributes of highlight, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{foldChange} only used in compareClusterResult object, fold Change of nodes, the default value is NULL. 
##'         If the user provides the Fold Change value of the nodes, 
##'         it can be used to set the size of the gene node.
##'         \item \code{category_node} Number indicating the amount by which plotting category
##'         nodes should be scaled relative to the default, the default value is 1.
##'         \item \code{gene_node} Number indicating the amount by which plotting gene nodes
##'         should be scaled relative to the default, the default value is 1.
##'         \item \code{category_label} Scale of category node label size, the 
##'         default value is 1.
##'         \item \code{gene_label} Scale of gene node label size, the default value is 1.
##'     }
##' @param hilight.params list, the parameters to control the attributes of highlighted nodes and edges.
##' see the hilight.params in the following.
##' hilight.params control the attributes of highlight, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{category} category nodes to be highlight.
##'         \item \code{alpha_hilight} alpha of highlighted nodes.
##'         \item \code{alpha_no_hilight} alpha of unhighlighted nodes.
##'     }
##' @importFrom ggraph geom_edge_arc
##' @importFrom ggplot2 scale_colour_gradient2
##' @importFrom rlang enquo
##' @author Guangchuang Yu
cnetplot.enrichResult <- function(x,
                     showCategory = 5,
                     foldChange = NULL,
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
                     color.params=list(
                         foldChange = NULL,
                         edge = FALSE,           
                         category = "#E5C494",   
                         gene = "#B3B3B3"        
                    ),
                     cex.params=list(
                         category_node = 1,      
                         gene_node = 1,          
                         category_label = 1,     
                         gene_label = 1         
                    ),
                     hilight.params=list(        
                         category = NULL,      
                         alpha_hilight = 1,      
                         alpha_no_hilight = 0.3  
                     ),
                     ...) {

    label_size_category <- 5
    label_size_gene <- 5
    node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
    
    params_df <- as.data.frame(rbind(
        c("foldChange", "color.params", "foldChange"),
        c("colorEdge", "color.params", "edge"),
        c("color_category", "color.params", "category"),
        c("color_gene", "color.params", "gene"),

        c("cex_category", "cex.params", "category_node"),
        c("cex_gene", "cex.params", "gene_node"),
        c("cex_label_category", "cex.params", "category_label"),
        c("cex_label_gene", "cex.params", "gene_label"))
    )
    colnames(params_df) <- c("original", "listname", "present")
    rownames(params_df) <- params_df$original


    
    default.color.params <- list(
        foldChange = NULL,
        edge = FALSE,           
        category = "#E5C494",   
        gene = "#B3B3B3"        
    )
    default.cex.params <- list(
        category_node = 1,       
        gene_node = 1,                
        category_label = 1,      
        gene_label = 1        
    )
    default.hilight.params <- list(
        category = NULL,
        alpha_hilight = 1,
        alpha_no_hilight = 0.3
    )

    # use modifyList to change the values of parameter 
    color.params <- modifyList(default.color.params, color.params)
    cex.params <- modifyList(default.cex.params, cex.params)
    hilight.params <- modifyList(default.hilight.params, hilight.params)
    params_list <- list(x = x,
            showCategory = showCategory,
            foldChange = foldChange,
            layout = layout,
            colorEdge = colorEdge,
            circular = circular,
            node_label = node_label,
            cex_category = cex_category,
            cex_gene = cex_gene,
            cex_label_category = cex_label_category,
            cex_label_gene = cex_label_gene,
            color_category = color_category,
            color_gene = color_gene,
            shadowtext = shadowtext,
            color.params = color.params, 
            cex.params = cex.params,
            hilight.params = hilight.params
    )

    # get all parameters value
    args <- as.list(match.call())
    removed_params <- intersect(params_df$original, names(args))
    if (length(removed_params) > 0) {
        for (i in removed_params) {
            params_list[[params_df[i, 2]]][[params_df[i, 3]]] <- get(i)
            warn <- get_param_change_message(i, params_df)
            warning(warn)
        }
    }

    
    color.params <- params_list[["color.params"]]
    cex.params <- params_list[["cex.params"]]
    hilight.params <- params_list[["hilight.params"]]
 
    foldChange <- color.params[["foldChange"]]                               
    colorEdge <- color.params[["edge"]]
    color_category <- color.params[["category"]]
    color_gene <- color.params[["gene"]]


    cex_category <- cex.params[["category_node"]]
    cex_gene <- cex.params[["gene_node"]]
    cex_label_category <- cex.params[["category_label"]]
    cex_label_gene <- cex.params[["gene_label"]]


    hilight_category <- hilight.params[["category"]]
    alpha_hilight <- hilight.params[["alpha_hilight"]]
    alpha_nohilight <- hilight.params[["alpha_no_hilight"]]

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

    # add edge alpha
    hilight_category <- intersect(hilight_category, names(geneSets))

    if (!is.null(hilight_category) && length(hilight_category) > 0) {
        edges <- attr(E(g), "vnames")
        E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
        hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
        hilight_gene <- edges[hilight_edge]
        hilight_gene <- gsub(".*\\|", "", hilight_gene)
        E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
    } else {
        E(g)$alpha <- rep(0.8, length(E(g)))
    }

    show_legend <- c(FALSE, TRUE)
    names(show_legend) <- c("alpha", "color") 
    if (colorEdge) {
        E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
        edge_layer <- geom_edge(aes_(color = ~category, alpha = ~I(alpha)), 
            show.legend = show_legend)
    } else {
        edge_layer <- geom_edge(aes_(alpha = ~I(alpha)), colour='darkgrey',
            show.legend = FALSE)
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
 
        p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
       
        alpha_category <- c(rep(1, n), rep(0, nrow(p$data)-n))
        alpha_gene <- c(rep(0, n), rep(1, nrow(p$data)-n))

        if (!is.null(hilight_category) && length(hilight_category) > 0) {
            alpha_category <- c(rep(alpha_nohilight, n), rep(0, nrow(p$data)-n))
            alpha_gene <- c(rep(0, n), rep(alpha_nohilight, nrow(p$data)-n))
            alpha_gene[match(hilight_gene, p$data$name)] <- alpha_hilight
            alpha_gene[match(hilight_category, p$data$name)] <- alpha_hilight
        }

        p <- p + edge_layer +
            geom_node_point(aes_(size=~size), color=I(color_category),
                        data = NULL, show.legend = show_legend,
                        alpha = I(alpha_category)) +
            ggnewscale::new_scale_color() +
            geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size),
                data = NULL, alpha = I(alpha_gene)) +
            scale_size(range=c(3, 8) * cex_category) +  
            # scale_colour_gradient2(name = "fold change") + 
            set_enrichplot_color(colors = get_enrichplot_color(3), name = "fold change")
            

    } else {
        V(g)$color <- color_gene
        V(g)$color[1:n] <- color_category
        p <- ggraph(g, layout=layout, circular=circular)
        p$data[-(1:n), "size"] <- 3 * cex_gene
        p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
        p <- p + edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size, alpha=~I(alpha)))+
            scale_size(range=c(3, 8) * cex_category) 

    }

    p <- p + theme_void()

    if (node_label == "category") { 
        p$data[-c(1:n), "name"] <- NA          
        p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category)
    } else if (node_label == "gene") {
        p$data[1:n, "name"] <- NA 
        p <- add_node_label(p = p, data = NULL, label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
    } else if (node_label == "all") {
        p <- add_node_label(p = p, data = NULL,
            label_size_node = c(rep(label_size_category, n), rep(label_size_gene, nrow(p$data)-n)),
            cex_label_node = c(rep(cex_label_category, n), rep(cex_label_gene, nrow(p$data)-n)), 
            shadowtext = shadowtext_gene)
    }
    if (!is.null(foldChange)) {
        p <- p + guides(size  = guide_legend(order = 1), 
                        color = guide_colorbar(order = 2))
    }
    return(p + guides(alpha = "none"))
}


##' @param split Separate result by 'category' variable.
##' @param pie Proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'.
##' Will be removed in the next version.
##' @param legend_n Number of circle in legend, the default value is 5.
##' Will be removed in the next version.
##' @param x_loc,y_loc The location of scatterpie legend.
##' Will be removed in the next version.
##' @param pie.params list, the parameters to control the attributes of pie nodes.
##' see the pie.params in the following.
##' pie.params control the attributes of pie nodes, it can be referred to the following parameters:
##'     \itemize{
##'         \item \code{pie} proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'.
##'         \item \code{legend_n} number of circle in legend.
##'         \item \code{legend_loc_x, legend_loc_y} The location of scatterpie legend.
##'     }
##' @importFrom ggraph geom_edge_arc
##' @noRd
cnetplot.compareClusterResult <- function(x,
                     showCategory = 5,
                     foldChange   = NULL,
                     layout = "kk",
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
                     pie.params = list(
                         pie = "equal",            
                         legend_n = 5,             
                         legend_loc_x = NULL,      
                         legend_loc_y = NULL       
                     ),
                     cex.params=list(
                         foldChange = NULL,
                         category_node = 1,        
                         gene_node = 1,            
                         category_label = 1,       
                         gene_label = 1           
                    ),
                     hilight.params=list(
                         category = NULL,
                         alpha_hilight = 1,
                         alpha_no_hilight = 0.3
                     ),  
                     ...) {

    label_size_category <- 2.5
    label_size_gene <- 2.5
    range_category_size <- c(3, 8)
    range_gene_size <- c(3, 3)

    # change parameter name
    ##############################################################
    params_df <- as.data.frame(rbind(
        c("pie", "pie.params", "pie"),
        c("legend_n", "pie.params", "legend_n"),
        c("x_loc", "pie.params", "legend_loc_x"),
        c("y_loc", "pie.params", "legend_loc_y"),

        c("foldChange", "cex.params", "foldChange"),
        c("cex_category", "cex.params", "category_node"),
        c("cex_gene", "cex.params", "gene_node"),
        c("cex_label_category", "cex.params", "category_label"),
        c("cex_label_gene", "cex.params", "gene_label"))
    )
    colnames(params_df) <- c("original", "listname", "present")
    rownames(params_df) <- params_df$original

    default.pie.params <- list(
        pie = "equal",          
        legend_n = 5,             
        legend_loc_x = NULL,      
        legend_loc_y = NULL           
    )
    default.cex.params <- list(
        foldChange = NULL,     
        category_node = 1,      
        gene_node = 1,          
        category_label = 1,     
        gene_label = 1               
    )
    default.hilight.params <- list(
        category = NULL,
        alpha_hilight = 1,
        alpha_no_hilight = 0.3
    )

    # use modifyList to change the values of parameter 
    pie.params <- modifyList(default.pie.params, pie.params)
    cex.params <- modifyList(default.cex.params, cex.params)
    hilight.params <- modifyList(default.hilight.params, hilight.params)

    params_list <- list(x = x,
        showCategory = showCategory,
        foldChange = foldChange,
        layout = layout,
        circular = circular,
        node_label = node_label ,
        split = split,
        pie = pie,
        cex_category = cex_category,
        cex_gene = cex_gene,
        legend_n = legend_n,
        x_loc = x_loc,
        y_loc = y_loc,
        cex_label_category = cex_label_category,
        cex_label_gene = cex_label_gene,
        shadowtext = shadowtext,
        pie.params = pie.params,
        cex.params = cex.params,
        hilight.params = hilight.params  
    )

    # get all parameters value
    args <- as.list(match.call())
    removed_params <- intersect(params_df$original, names(args))
    if (length(removed_params) > 0) {
        for (i in removed_params) {
            params_list[[params_df[i, 2]]][[params_df[i, 3]]] <- get(i)
            warn <- get_param_change_message(i, params_df)
            warning(warn)
        }
    }

    pie.params <- params_list[["pie.params"]]
    cex.params <- params_list[["cex.params"]]
    hilight.params <- params_list[["hilight.params"]]

    pie <- pie.params[["pie"]]
    legend_n <- pie.params[["legend_n"]]
    x_loc <- pie.params[["legend_loc_x"]]
    y_loc <- pie.params[["legend_loc_y"]]


    foldChange <- cex.params[["foldChange"]]
    cex_category <- cex.params[["category_node"]]
    cex_gene <- cex.params[["gene_node"]]
    cex_label_category <- cex.params[["category_label"]]
    cex_label_gene <- cex.params[["gene_label"]]


    # default.hilight.params <- list(
    #     category = NULL,
    #     alpha_hilight = 1,
    #     alpha_no_hilight = 0.3
    # )
    # hilight.params <- reset_params(defaultp=default.hilight.params, 
    #                                inputp=enquo(hilight.params))
    hilight_category <- hilight.params[["category"]]
    alpha_hilight <- hilight.params[["alpha_hilight"]]
    alpha_nohilight <- hilight.params[["alpha_no_hilight"]]

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
    hilight_category <- intersect(hilight_category, y$Description)
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
    edge_layer <- geom_edge(aes_(alpha = ~I(alpha)), colour='darkgrey',
            show.legend = FALSE)
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
 
        p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
        p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
        p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category)   
        return(p)
    }
 
    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        p <- ggraph(g, "tree") + edge_layer
    } else if (length(unique(y$Cluster)) == 1) {
        color_category = "#E5C494"
        color_gene = "#B3B3B3"
        size <- sapply(geneSets, length)
        V(g)$size <- min(size)/2
        n <- length(geneSets)
        V(g)$size[1:n] <- size
        V(g)$color <- color_gene
        V(g)$color[1:n] <- color_category

        if (!is.null(hilight_category) && length(hilight_category) > 0) {
            edges <- attr(E(g), "vnames")
            E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
            hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
            hilight_gene <- edges[hilight_edge]
            hilight_gene <- gsub(".*\\|", "", hilight_gene)
            E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
        } else {
            E(g)$alpha <- rep(0.8, length(E(g)))
        }
        p <- ggraph(g, layout=layout, circular=circular)
        p$data[-(1:n), "size"] <- 3 * cex_gene
        p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
        p <- p + edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size, alpha=~I(alpha)))+
            scale_size(range=c(3, 8) * cex_category) 
        
    } else {
        if (!is.null(hilight_category) && length(hilight_category) > 0) {
            edges <- attr(E(g), "vnames")
            E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
            hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
            hilight_gene <- edges[hilight_edge]
            hilight_gene <- gsub(".*\\|", "", hilight_gene)
            E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
        } else {
            E(g)$alpha <- rep(0.8, length(E(g)))
        }
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

            p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
            ID_Cluster_mat2$alpha <-  p$data$alpha       
            p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius,alpha=~I(alpha)),
                    data=ID_Cluster_mat2[1:n, ],
                    cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-4)], 
                    color=NA) +
                geom_scatterpie_legend(ID_Cluster_mat2$radius[1:n],
                    x=x_loc, y=y_loc + 3, n = legend_n, labeller=function(x) round(x^2 * sum_yunion / cex_category))  +
                geom_scatterpie(aes_(x=~x,y=~y,r=~radius,alpha=~I(alpha)),
                    data=ID_Cluster_mat2[-(1:n), ],
                    cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-4)],
                    color=NA,
                    show.legend = FALSE) +
                coord_equal()+
                geom_scatterpie_legend(ID_Cluster_mat2$radius[(n+1):nrow(ID_Cluster_mat2)],
                    x=x_loc, y=y_loc, n = legend_n,
                    labeller=function(x) round(x*2/(min(sizee))/sqrt(cex_gene),1)) +
                ggplot2::annotate("text", x = x_loc + 3, y = y_loc, label = "log2FC")  +
                ggplot2::annotate("text", x = x_loc + 3, y = y_loc + 3, label = "gene number")

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
        p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight) 
        ID_Cluster_mat2$alpha <-  p$data$alpha  
        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius,alpha=~I(alpha)),
                data=ID_Cluster_mat2[1:n, ],
                cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-4)], color=NA) +
            geom_scatterpie(aes_(x=~x,y=~y,r=~radius,alpha=~I(alpha)),
                data=ID_Cluster_mat2[-(1:n), ],
                cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-4)],
                color=NA, show.legend = FALSE) +
            coord_equal() +
            geom_scatterpie_legend(ID_Cluster_mat2$radius[1:n],
                    x=x_loc, y=y_loc, n = legend_n, labeller=function(x) round(x^2 * sum_yunion / cex_category)) +
            ggplot2::annotate("text", x = x_loc + 3, y = y_loc, label = "gene number")
        ## add node label
        p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
            cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
        p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext_category)
        p <- p + theme_void() + labs(fill = "Cluster")
        return(p)
    }
    title <- colnames(ID_Cluster_mat2)[1]
    # V(g)$size <- ID_Cluster_mat2$radius
    # V(g)$color <- "#B3B3B3"
    # V(g)$color[1:n] <- "#E5C494"  
    

    p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
    p <- add_node_label(p = p, data = p$data[-c(1:n),], label_size_node = label_size_gene,
        cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
    p <- add_node_label(p = p, data = p$data[1:n,], label_size_node = label_size_category,
        cex_label_node = cex_label_category, shadowtext = shadowtext_category) + theme_void()
    if (length(unique(y$Cluster)) > 1) {
        p <- p + theme(legend.position="none")
    }
    p
}