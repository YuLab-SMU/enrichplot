##' Functional grouping network diagram for enrichment result of
##' over-representation test or gene set enrichment analysis
##' @importFrom igraph layout_with_fr
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_discrete
##' @importFrom ggplot2 scale_size_continuous
##' @importFrom stats kmeans
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @importFrom GOSemSim godata
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue, p.adjust or qvalue
##' @param with_edge if TRUE, draw the edges of the network diagram
##' @param line_scale scale of line width
##' @param method method of enrichment, one of "enrichGO", "enrichDO" and "other" 
##' @param N the number of words in the cluster tags 
##' @param ncluster the number of clusters 
##' @param OrgDb OrgDb object 
##' @param split separate result by 'category' variable
##' @param ont one of 'BP', 'MF', 'CC'
##' @param measure method for calculating semantic similarity, one of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods
##' @param min_edge minimum percentage of overlap genes to display the edge, should between 0 and 1, default value is 0.2
##' @export
##'
##' @examples
##' \dontrun{
##'     library(clusterProfiler)
##'     data(geneList, package = "DOSE")
##'	    gene <- names(geneList)[abs(geneList) > 2]
##'	    ego <- enrichGO(gene  = gene,
##'         universe      = names(geneList),
##'         OrgDb         = org.Hs.eg.db,
##'         ont           = "CC",
##'         pAdjustMethod = "BH",
##'         pvalueCutoff  = 0.01,
##'         qvalueCutoff  = 0.05,
##'         readable      = TRUE)
##'	    emapplot_cluster(ego)
##'    }
emapplot_cluster <- function(x, showCategory = nrow(x), color = "p.adjust", line_scale = 0.1, with_edge = FALSE,
     method = "other", N = 4, ncluster = NULL, OrgDb = NULL, split=NULL, ont="BP", measure = "Wang", min_edge=0.2){
    
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
    
    
    g <- get_igraph(x=x, y=y, n=n, color=color, line_scale=line_scale, min_edge=min_edge)
    if(n == 1) {
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    }
    edgee <- igraph::get.edgelist(g)
    ## Calculate the semantic similarity or overlap between two nodes    
    y_data <- as.data.frame(y)
    rownames(y_data) <- y_data$Description
    edgee2 <- data.frame(GO1 = y_data[edgee[, 1], "ID"], 
        GO2 = y_data[edgee[, 2], "ID"], stringsAsFactors = FALSE)
    

    if(method == "enrichGO") {
        d <- godata(OrgDb = OrgDb, ont=ont, computeIC=FALSE)
        fun1 <- function(term_id) {
            GOSemSim::goSim(term_id[1], term_id[2], semData=d, measure=measure)
        } 
    } else if(method == "enrichDO") {
    
        fun1 <- function(term_id) {
            DOSE::doSim(term_id[1], term_id[2], measure=measure)
        } 
    } else {
        fun1 <- function(term_id) {   
            id1 <- y[term_id[1], "geneID"]
            id1_1 <- unlist(strsplit(id1, split="/"))
            id2 <- y[term_id[2], "geneID"]
            id2_2 <- unlist(strsplit(id2, split="/"))
            overlap_ratio(id1_1, id2_2)
        }
    }
         
    ## Set layout according to semantic similarity or overlap
    edge_w <- apply(edgee2, 1, fun1)
    set.seed(123)   
    lw <- layout_with_fr(g, weights=edge_w)
    
    p <- ggraph::ggraph(g, layout=lw)
    # cluster_label1 <- lapply(clusters, function(i){i[order(y[i, "pvalue"])[1]]})
    
    ## Using k-means clustering to group
    pdata2 <- p$data
    dat <- data.frame(x=pdata2$x, y=pdata2$y)
    if(is.null(ncluster)){
        pdata2$color <- kmeans(dat, ceiling(sqrt(nrow(dat))))$cluster
    } else {
        if(ncluster > nrow(dat)) ncluster <- nrow(dat)
        pdata2$color <- kmeans(dat, ncluster)$cluster
    }
    
    goid <- y_data$ID
    cluster_color <- unique(pdata2$color)
    clusters <- lapply(cluster_color, function(i){goid[which(pdata2$color == i)]})
    cluster_label <- sapply(cluster_color,  wordcloud_i, pdata2 = pdata2, N=N)
    names(cluster_label) <- cluster_color
    pdata2$color <- cluster_label[as.character(pdata2$color)]
    p$data <- pdata2
    p <- p + ggforce::geom_mark_hull(aes_(label=~color, x=~x, y=~y, colour=~color), show.legend = FALSE) +
        geom_node_point(aes_(colour=~color, size=~size)) + scale_color_discrete(name="cluster")  + 
        scale_size_continuous(name = "number of genes")
    if(with_edge) {
        p <- p + ggraph::geom_edge_link(alpha=.8, aes_(width=~I(width*line_scale)), colour='darkgrey')
    }
    return(p)

}
