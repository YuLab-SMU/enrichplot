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
##'                   one of 'category', 'gene', 'all' and 'none', default is "all".
##' @importFrom ggraph geom_edge_arc
##' @author Guangchuang Yu
cnetplot.enrichResult <- function(x,
                     showCategory = 5,
                     foldChange   = NULL,
                     layout = "kk",
                     colorEdge = FALSE,
                     circular = FALSE,
                     node_label = "all",
                     ...) {

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
        palette <- fc_palette(fc)
        p <- ggraph(g, layout=layout, circular = circular) +
            edge_layer +
            geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
            scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494")
    } else {
        V(g)$color <- "#B3B3B3"
        V(g)$color[1:n] <- "#E5C494"
        p <- ggraph(g, layout=layout, circular=circular) +
            edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size))
    }

    p <- p + scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
        theme_void()


    if (node_label == "category") {
        p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,])
    } else if (node_label == "gene") {
        p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),], repel=TRUE)
    } else if (node_label == "all") {
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE)
    } 

    return(p)
}

##' @param colorEdge whether coloring edge by enriched terms
##' @param circular whether using circular layout
##' @param node_label select which labels to be displayed.
##'                   one of 'category', 'gene', 'all' and 'none', default is "all".
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @param pie_scale scale of pie plot
##' @param legend_n number of circle in legend
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
                     pie_scale = 1,
                     legend_n = 5,
                     ...) {
    n <- update_n(x, showCategory)
    if (n == 0) {
        stop("no enriched term found...")
    }
  
    y <- fortify.compareClusterResult(x, showCategory=showCategory,
                                      includeAll=TRUE, split=split)
    y$Cluster = sub("\n.*", "", y$Cluster)
    
   
    #n <- update_n(x, showCategory)  
    y_union <- merge_compareClusterResult(y)
    node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
    ## when y just have one line
   
    
    if (circular) {
        layout <- "linear"
        geom_edge <- geom_edge_arc
    } else {
        geom_edge <- geom_edge_link
    }

 
    #geneSets <- extract_geneSets(x, showCategory)
    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/", fixed = TRUE), y_union$ID)
    n <- length(geneSets)
    g <- list2graph(geneSets)
    edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
    if(is.null(dim(y)) | nrow(y) == 1) {
        V(g)$size <- 1
        V(g)$size[1] <- 3
        V(g)$color <- "#B3B3B3"
        V(g)$color[1] <- "#E5C494"
        p <- ggraph(g, layout=layout, circular=circular) +
            edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size))
        return(p)
    }
    
    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
        p <- ggraph(g) + edge_layer   
    } else {
        p <- ggraph(g, layout=layout, circular=circular) + edge_layer
    }

    
    #pie chart begin
    #obtain the cluster distribution of each GO term and gene
    ID_Cluster_mat <- deal_data_pie(y, pie=pie)

    gene_Cluster_mat <- deal_gene_pie(y)
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
    ii <- match(rownames(ID_Cluster_mat2)[1:n], y_union$ID)
    sizee <- sqrt(y_union[ii,9] / sum(y_union[ii,9])) * pie_scale
    ID_Cluster_mat2$radius <- min(sizee)/2
    ID_Cluster_mat2$radius[1:n] <- sizee
    x_loc1 <- min(ID_Cluster_mat2$x)
    y_loc1 <- min(ID_Cluster_mat2$y)
    #node_label
    if (node_label == "category") {
        p$data$name[(n+1):nrow(p$data)] <- ""
    } else if (node_label == "gene") {
        p$data$name[1:n] <- ""
    } else if (node_label == "none") {
        p$data$name <- ""
    }
    if(ncol(ID_Cluster_mat2) > 4) {
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
            ID_Cluster_mat2$radius <- min(sizee)/2*gene_fc2
            ID_Cluster_mat2$radius[1:n] <- sizee
            p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat2,
                                 cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)],color=NA) +
                coord_equal()+
                geom_scatterpie_legend(ID_Cluster_mat2$radius[(n+1):nrow(ID_Cluster_mat2)], x=x_loc1, y=y_loc1, 
                n = legend_n, labeller=function(x) round(x*2/(min(sizee)),3)) +
                geom_node_text(aes_(label=~name), repel=TRUE, size=2.5) + theme_void() +
                labs(fill = "Cluster")
            return(p)        
        }
        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat2,
                                 cols=colnames(ID_Cluster_mat2)[1:(ncol(ID_Cluster_mat2)-3)],color=NA) +
        coord_equal()+
        geom_node_text(aes_(label=~name), repel=TRUE, size=2.5) + theme_void() +
        labs(fill = "Cluster")
        return(p)    
    }
    title <- colnames(ID_Cluster_mat2)[1]
    V(g)$size <- ID_Cluster_mat2$radius
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"

    ggraph(g, layout=layout, circular=circular) + 
    edge_layer + geom_node_point(aes_(color=~I(color), size=~size)) + labs(title= title) +
    scale_size(range=c(3, 8) * pie_scale)
}

##' Prepare the data for the pie plot
##'
##' @param y a data.frame of clusterProfiler result
##' @return a data.frame
##' @noRd
deal_gene_pie <- function(y) {

    #data_pie <- as.matrix(y[,c(1,2,10)])
    gene_pie <- as.matrix(y[,c(1,2,9)])
    gene_pie2 <- NULL
    for(i in seq_len(dim(gene_pie)[1])){
        gene <- unlist(strsplit(gene_pie[i, 3], "/"))
        for(j in seq_len(length(gene))) {
            
            gene_pie2 <- rbind(gene_pie2,c(gene[j],gene_pie[i,1:2]))
        }
    }
    gene_pie2 <- unique(gene_pie2)
    genes_unique <- unique(gene_pie2[,1])
    Cluster_unique <- unique(gene_pie2[,2])

    #ID_Cluster_mat <- matrix(0,length(ID_unique),length(Cluster_unique))
    gene_Cluster_mat <- matrix(0,length(genes_unique),length(Cluster_unique))   
    #rownames(ID_Cluster_mat) <- ID_unique
    #colnames(ID_Cluster_mat) <- Cluster_unique
    rownames(gene_Cluster_mat) <- genes_unique
    colnames(gene_Cluster_mat) <- Cluster_unique
    
    #ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
    gene_Cluster_mat <- as.data.frame(gene_Cluster_mat, stringAsFactors = FALSE)
    for(i in seq_len(nrow(gene_pie2))) {
        gene_Cluster_mat[gene_pie2[i,1],gene_pie2[i,2]] <- 1
    }
    return(gene_Cluster_mat)
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
    ldf <- lapply(1:length(inputList), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}


