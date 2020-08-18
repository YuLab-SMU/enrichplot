##' Functional grouping network diagram for enrichment result of
##' over-representation test or gene set enrichment analysis
##' @importFrom igraph layout_with_fr
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_discrete
##' @importFrom ggplot2 scale_size_continuous
##' @importFrom ggplot2 scale_fill_discrete
##' @importFrom stats kmeans
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @importFrom GOSemSim godata
##' @importFrom shadowtext geom_shadowtext
##' @importFrom ggnewscale new_scale_fill
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue, p.adjust or qvalue
##' @param with_edge if TRUE, draw the edges of the network diagram
##' @param line_scale scale of line width
##' @param method method of calculating the similarity between nodes, one of "goSim", "doSim" and "other" 
##' @param nWords the number of words in the cluster tags 
##' @param nCluster the number of clusters 
##' @param OrgDb OrgDb object 
##' @param split separate result by 'category' variable
##' @param ont one of 'BP', 'MF', 'CC'
##' @param measure method for calculating semantic similarity, one of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods
##' @param min_edge minimum percentage of overlap genes to display the edge, should between 0 and 1, default value is 0.2
##' @param cluster_label_scale scale of cluster labels size
##' @export
##'
##' @examples
##' \dontrun{
##'     library(clusterProfiler)
##'     data(geneList, package = "DOSE")
##'        gene <- names(geneList)[abs(geneList) > 2]
##'        ego <- enrichGO(gene  = gene,
##'         universe      = names(geneList),
##'         OrgDb         = org.Hs.eg.db,
##'         ont           = "CC",
##'         pAdjustMethod = "BH",
##'         pvalueCutoff  = 0.01,
##'         qvalueCutoff  = 0.05,
##'         readable      = TRUE)
##'      emapplot_cluster(ego)
##'    }
emapplot_cluster <- function(x, showCategory = nrow(x), color = "p.adjust", line_scale = 0.1, with_edge = TRUE,
     method = "other", nWords = 4, nCluster = NULL, OrgDb = NULL, split = NULL, ont = "BP", measure = "Wang", 
     min_edge = 0.2, cluster_label_scale = 1){
    
    n <- update_n(x, showCategory)
    if(class(x) == "compareClusterResult") {
        y <- fortify(x, showCategory=showCategory,
                                      includeAll=TRUE, split=split)
        y$Cluster = sub("\n.*", "", y$Cluster)
    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    ## Data structure transformation, combining the same ID (Description) genes       
        y <- merge_compareClusterResult(y)        
    } else {
        y <- as.data.frame(x)
    }
    
    
    g <- get_igraph(x=x, y=y, n=n, color=color, line_scale=line_scale, min_edge=min_edge, 
        method = method, measure = measure, OrgDb = OrgDb, ont = ont)
    if(n == 1) {
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    }
    edgee <- igraph::get.edgelist(g)
    ## Calculate the semantic similarity or overlap between two nodes    
    y_data <- as.data.frame(y)
    rownames(y_data) <- y_data$Description
    edgee2 <- data.frame(GO1 = y_data[edgee[, 1], "ID"], 
        GO2 = y_data[edgee[, 2], "ID"], stringsAsFactors = FALSE)
    
    edge_w <- E(g)$weight
    set.seed(123)   
    lw <- layout_with_fr(g, weights=edge_w)
    
    p <- ggraph::ggraph(g, layout=lw)
    # cluster_label1 <- lapply(clusters, function(i){i[order(y[i, "pvalue"])[1]]})
    
    ## Using k-means clustering to group
    pdata2 <- p$data
    dat <- data.frame(x = pdata2$x, y = pdata2$y)
    colnames(pdata2)[5] <- "color2"
    
    if(is.null(nCluster)){
        pdata2$color <- kmeans(dat, ceiling(sqrt(nrow(dat))))$cluster
    } else {
        if(nCluster > nrow(dat)) nCluster <- nrow(dat)
        pdata2$color <- kmeans(dat, nCluster)$cluster
    }
    
    goid <- y_data$ID
    cluster_color <- unique(pdata2$color)
    clusters <- lapply(cluster_color, function(i){goid[which(pdata2$color == i)]})
    cluster_label <- sapply(cluster_color,  wordcloud_i, pdata2 = pdata2, nWords=nWords)
    names(cluster_label) <- cluster_color
    pdata2$color <- cluster_label[as.character(pdata2$color)]
    p$data <- pdata2
    ## Take the location of each group's center nodes as the location of the label
    label_x <- stats::aggregate(x ~ color, pdata2, mean)
    label_y <- stats::aggregate(y ~ color, pdata2, mean)
    label_location <- data.frame(x = label_x$x, y = label_y$y, label = label_x$color)
    
    show_legend <- c(TRUE, FALSE)
    names(show_legend) <- c("fill", "color")
    
    if(with_edge) {
        p <-  p +  ggraph::geom_edge_link(alpha = .8, aes_(width =~ I(width*line_scale)), colour='darkgreen')
    }
    
    
    p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color, fill =~ color), show.legend = show_legend) + 
        scale_fill_discrete(name = "cluster") + new_scale_fill() + 
        geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ color2, size =~ size)) + 
        scale_size_continuous(name = "number of genes")  + 
        scale_fill_continuous(low = "red", high = "blue", name = color, guide = guide_colorbar(reverse = TRUE))  + 
        geom_shadowtext(data = label_location, aes_(x =~ x, y =~ y, label =~ label), 
            size = 5 * cluster_label_scale, check_overlap = TRUE)
        
}
