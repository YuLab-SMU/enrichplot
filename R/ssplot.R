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
##' @importFrom ggplot2 theme_classic
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @param color Variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @param cex_line Scale of line width
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##' @param cex_label_category Scale of category node label size.
##' @param cex_category Number indicating the amount by which plotting category
##' nodes should be scaled relative to the default.
##' @param shadowtext a logical value, whether to use shadow font.
##' @param group_category a logical, if TRUE(the default), group the category.
##' @param node_label Select which labels to be displayed,
##' one of 'category', 'group', 'all' and 'none'.
##' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
##' The default is FALSE.
##' @param label_style style of group label, one of "shadowtext" and "ggforce".
##' @param repel whether to correct the position of the label. Defaults to FALSE.
##' @param cex_label_group Numeric, scale of group labels size, the default value is 1.
##' @param nWords Numeric, the number of words in the cluster tags, the default value is 4.
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' @param method The function used for MDS, 
##' one of "cmdscale" (the default), "sammon", "monoMDS" and "isoMDS".
##' 
##' additional parameters can refer the following parameters.
##'     \itemize{
##'        \item \code{force} Force of repulsion between overlapping text labels. Defaults to 1. 
##'        \item \code{nudge_x, nudge_y} Horizontal and vertical adjustments to nudge 
##'         the starting position of each text label. 
##'        \item \code{direction} "both", "x", or "y" – direction in which to adjust position of labels.
##'     }
##' 
ssplot.enrichResult <- function(x, showCategory = 30, color="p.adjust",
                                min_edge=0.2,
                                cex_label_category  = 1, cex_category = 1,
                                cex_line = 1, shadowtext = TRUE,
                                group_category = FALSE,
                                node_label  = "category",
                                with_edge = FALSE, label_format = 30,
                                group_legend = FALSE, 
                                label_style = "shadowtext", repel = FALSE,
                                cex_label_group = 1, nWords = 4, 
                                nCluster = NULL, method = "cmdscale", ...) {
    has_pairsim(x)
    label_size_category <- 5
    label_group <- 3
    n <- update_n(x, showCategory)
    y <- as.data.frame(x)
    ## get graph.data.frame object
    g <- get_igraph(x=x, y=y, n=n, color=color, cex_line=cex_line,
                    min_edge=min_edge)
    ## If there is only one point, then add a dot and label, then return directly.
    if(n == 1) {
        return(ggraph(g,"tree") + geom_node_point(color="red", size=5) +
               geom_node_text(aes_(label=~name)))
    }
    ## get ggraph object
    p <- ggraph(g, layout= "nicely")

    ## Then the similarity matrix is used to reduce the 
    ## dimension and redistribute the node coordinates.
    ## Complete similarity matrix
    if (is.numeric(n)) {
        keep <- seq_len(n)
    } else {
        keep <- match(n, rownames(x@termsim))
    }
    if (length(keep) == 0) {
        stop("no enriched term found...")
    }
    termsim2 <- fill_termsim(x, keep)

    ## Dimensionality reduction and K-means clustering using similarity matrix
    data <- get_cluster(termsim2, nCluster = nCluster, method = method)

    ## Add other information of the enrichment result to the data object.
    terms <- rownames(termsim2)
    object <- as.data.frame(x)
    rownames(object) <- object$Description
    data <- cbind(data, object[terms, ])
    
    ## Reset node coordinates according to the results of clustering
    pdata2 <- p$data
    pdata2$x <- data$x
    pdata2$y <- data$y

    ## Set axis label according to method
    p <- adj_axis(p = p, data = data)

    ## The next order: 
    ## (1) add edges
    ## (2) add circles, 
    ## (3) add dots
    ## (4) add dots labels
    ## (5 )add circles labels
    ## add edge
    p$data <- pdata2  
    if (with_edge & length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }
    

    # if show group cricle or group label, Process p$data and assign color to the group label
    if (group_category || node_label == "all" || node_label == "group") {     
        pdata2$color2 <- data$color
        pdata2 <- groupNode(pdata2, y = y, nWords = nWords, nCluster = nCluster)
        p$data <- pdata2
    }

    ## if group_category, add circles
    if (group_category) {
        p <- add_ellipse(p = p, group_legend = group_legend, 
            label_style = label_style)
    }


    ## add dot
    p <- add_category_nodes(p = p, cex_category = cex_category, color)  
    if (node_label == "all" || node_label == "category") 
        p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
            cex_label_node = cex_label_category, shadowtext = shadowtext)

    if (node_label == "all" || node_label == "group") {   
        label_location <- get_label_location(pdata2, label_format)
        p <- add_group_label(repel = repel, shadowtext = shadowtext, p = p,
            label_location = label_location, label_group = label_group,
            cex_label_group = cex_label_group)
    }
    return(p + theme_classic() + coord_equal())
}




##' @rdname ssplot
##' @importFrom igraph E "E<-"
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_classic
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
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'
##' @param legend_n number of circle in legend
##' @importFrom stats setNames
ssplot.compareClusterResult <- function(x, showCategory = 30,
                                        split = NULL, pie = "equal",
                                        legend_n = 5, cex_category = 1,
                                        cex_line = 1, min_edge=0.2, 
                                        cex_label_category  = 1, 
                                        shadowtext = TRUE, 
                                        with_edge = FALSE,
                                        group_category = FALSE, 
                                        label_format = 30,
                                        group_legend = FALSE,
                                        node_label  = "category",
                                        label_style = "shadowtext", 
                                        repel = FALSE, cex_label_group = 1,
                                        nWords = 4, nCluster = NULL, 
                                        method = "cmdscale", ...) {
                                       
    has_pairsim(x)
    label_size_category <- 3
    label_group <- 3
    y <- get_selected_category(showCategory, x, split)
    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- merge_compareClusterResult(y)
     
    ## get ggraph object and add edge, here the color parameter is useless, 
    ## just to match ssplot.enrichResult
    p <- build_ggraph(x = x, y = y, y_union = y_union, cex_category = cex_category, 
        pie = pie, layout = "nicely", cex_line=cex_line,
                        min_edge=min_edge, pair_sim = x@termsim,
                        method = x@method, with_edge = with_edge)
    if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(y_union)) | nrow(y_union) == 1)
        return(p)
    ## Then the similarity matrix is used to reduce the 
    ## dimension and redistribute the node coordinates.
    ## Complete similarity matrix
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    termsim2 <- fill_termsim(x, rownames(ID_Cluster_mat))

    ## Dimensionality reduction and K-means clustering using similarity matrix
    data <- get_cluster(termsim2, nCluster = nCluster, method = method)

    ## Add other information of the enrichment result to the data object.
    data <- cbind(data, y_union)
    
    ## Reset node coordinates according to the results of clustering
    pdata2 <- p$data
    pdata2$x <- data$x
    pdata2$y <- data$y
    ## Set axis label according to the method parameter
    p <- adj_axis(p = p, data = data)
    ## The next order: add a circle, add a point, add a point label, add a circle label.
    # if show group cricle or group label, Process p$data and assign color to the group label  
    pdata2$color2 <- data$color  
    # if show group cricle or group label, Process p$data and assign color to the group label
    if (group_category || node_label == "all" || node_label == "group") {    
       pdata2 <- groupNode(pdata2, y = y, nWords = nWords, nCluster = nCluster)
    }   
    p$data <- pdata2   
    ## add circle
    if (group_category) {
        p <- add_ellipse(p = p, group_legend = group_legend, 
            label_style = label_style)
    }
    
    
    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- get_pie_data(y = y, pie = pie, y_union = y_union, 
        aes_axis = 0.0125, pdata2 = pdata2, cex_category = cex_category)

    ## get the location of legend for size of pie
    x_loc1 <- min(ID_Cluster_mat$x)
    y_loc1 <- min(ID_Cluster_mat$y)
    
    ## add dot and node label
    p <- add_pie_node(p = p, ID_Cluster_mat = ID_Cluster_mat, 
                  node_label = node_label, cex_category = cex_category,
                  aes_axis = 0.0125, 
                  cex_label_category = cex_label_category,
                  x_loc1 = x_loc1, y_loc1 = y_loc1,
                  shadowtext = shadowtext, legend_n = legend_n,
                  pdata2 = pdata2, 
                  label_size_category = label_size_category)  
    ## add group label
    if (node_label == "all" || node_label == "group") {   
        label_location <- get_label_location(pdata2, label_format)
        p <- add_group_label(repel = repel, shadowtext = shadowtext, p = p,
            label_location = label_location, label_group = label_group,
            cex_label_group = cex_label_group)
    }     
    return(p + theme_classic())
}


#' Group terms
#'
#' @param termsim2 Similarity matrix of terms.
#' @param nCluster Numeric, the number of cluster.
#'
#' @noRd
get_cluster <- function(termsim2 = termsim2, nCluster, method) {
    ## If the similarity between the two terms is 1, an error will be reported, so fine-tuning.
    termsim2[which(termsim2 == 1)] <- 0.99999
    for (i in seq_len(nrow(termsim2))) termsim2[i, i] <- 1
    
    dune.dist <- stats::as.dist(1- termsim2)
    # res <- ape::pcoa(dune.dist, correction="none")
    if (method == "cmdscale") {
        res <- stats::cmdscale(dune.dist, eig = T)   
    }

    if (method == "sammon") {
        y <- stats::cmdscale(dune.dist)
        ## If the matrix y has duplicate rows it will report an error, so perturb slightly  
        dup <- which(duplicated(y) == TRUE)  
        y[dup, 1] <- y[dup, 1] + 10^-7 * seq_len(length(dup))  
        res <- MASS::sammon(d = dune.dist, y = y)
    }  
    
    if (method == "monoMDS" || method == "isoMDS") {
        res <- vegan::metaMDS(dune.dist, method = method)
    }
    data <- as.data.frame(res$points[, 1:2])
    colnames(data) <- c('x', 'y')

    if (is.null(nCluster)){
        data$color <- stats::kmeans(data, floor(sqrt(nrow(data))))$cluster
    } else {
        if (nCluster > nrow(data)) nCluster <- nrow(data)
        data$color <- stats::kmeans(data, nCluster)$cluster
    }
    data$color <- as.character(data$color)
    if (method == "cmdscale") {
        data$pocas <- 0
        data$pocas[seq_len(length(res$eig))] <- as.numeric(res$eig)
    } else {
        data$stress <- res$stress
    }
    return(data)
}
