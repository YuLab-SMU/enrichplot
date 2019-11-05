##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "enrichResult"),
          function(x, showCategory = 30, color = "p.adjust", layout = "kk", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, color = "p.adjust", layout = "kk", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 5, color = "p.adjust", layout = "kk", ...) {

              emapplot.compareClusterResult(x, showCategory = showCategory, color=color,
                                            layout = layout, ...)
          })

emap_graph_build <- function(n,y,geneSets,color) {
    if (is.null(dim(y)) | nrow(y) == 1) {
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- as.character(y$Description)
        V(g)$color <- "red"
        ##return(ggraph(g))
        } else {
        id <- y[,"ID"]
        geneSets <- geneSets[id]
        n <- nrow(y) #
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$Description
        
        for (i in 1:n) {
            for (j in i:n) {
                w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
        
        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
        g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
        ## g <- delete.edges(g, E(g)[wd[,3] < 0.05])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
        
        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt
        
        colVar <- y[idx, color]
        V(g)$color <- colVar
    }
 
    return(g)
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
##' @author Guangchuang Yu
emapplot.enrichResult <- function(x, showCategory = 30, color="p.adjust", layout = "kk", ...) {
    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    if (is.numeric(n)) {
        y <- y[1:n,]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }


    if (n == 0) {
        stop("no enriched term found...")
    } 
    
    g <- emap_graph_build(n=n,y=y,geneSets=geneSets,color=color)
    if(n == 1) {
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))

    } else {
    p <- ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    } 
    p + geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8))
    }
}



merge_compareClusterResult <- function(yy) {
    yy_union<- yy[!duplicated(yy$ID),]
    yy_ids <- lapply(split(yy, yy$ID), function(x) {
        ids <- unique(unlist(strsplit(x$geneID, "/")))
        cnt <- length(ids)
        list(ID=paste0(ids, collapse="/"), cnt=cnt)
    })

    ids <- vapply(yy_ids, function(x) x$ID, character(1))
    cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))
    
    yy_union$geneID = ids[yy_union$ID]
    yy_union$Count = cnt[yy_union$ID]
    yy_union$Cluster = NULL
    yy_union
}

deal_data_pie <- function(y) {
    data_pie <- as.matrix(y[,1:2])
	ID_unique <- unique(data_pie[,2])
    Cluster_unique <- unique(data_pie[,1])
	ID_Cluster_mat <- matrix(0,length(ID_unique),length(Cluster_unique))
	rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    for(i in seq_len(nrow(data_pie))) {
        ID_Cluster_mat[data_pie[i,2],data_pie[i,1]] <- 1
    }
	return(ID_Cluster_mat)
}





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
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 coord_equal
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ylim
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @importFrom DOSE geneID
##' @importClassesFrom DOSE compareClusterResult
##' @param split separate result by 'category' variable
##' @method fortify compareClusterResult
##' @importFrom scatterpie geom_scatterpie
##' @importFrom stats setNames
##' @noRd
emapplot.compareClusterResult <- function(x, showCategory = 5, color = "p.adjust",
                                          layout = "kk", split=NULL, ...) {

    region <- radius <- NULL

    ## pretreatment of x, just like dotplot do
    y <- fortify.compareClusterResult(x, showCategory=showCategory, 
                                      includeAll=TRUE, split=split)

    y$Cluster = sub("\n.*", "", y$Cluster)


    ## geneSets <- geneInCategory(x) ## use core gene for gsea result

    ## Data structure transformation, combining the same ID (Description) genes
    n <- update_n(x, showCategory)
    
    y_union <- merge_compareClusterResult(y)
    
    if (n == 0) {
        stop("no enriched term found...")
    } 
    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/", fixed = TRUE), y_union$ID)
    g <- emap_graph_build(n=n,y=y_union,geneSets=geneSets,color=color)
    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
         ##return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
        p <- ggraph(g)
		ID_Cluster_mat <- deal_data_pie(y)
		ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1)
		colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],"x","y","radius")
		ID_Cluster_mat <- as.data.frame(ID_Cluster_mat)
        p + scatterpie::geom_scatterpie(aes(x=x,y=y,r=radius), data=ID_Cluster_mat,
                                    cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA)+
            xlim(-3,3) + ylim(-3,3) + coord_equal()+ geom_node_text(aes_(label=~name), repel=TRUE) + theme_void()
            									
    } else {
        p <- ggraph(g, layout=layout)
        if (length(E(g)$width) > 0) {
            p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
        }
     
        ## then add the pie plot    
        ## Get the matrix data for the pie plot
        ID_Cluster_mat <- deal_data_pie(y)
    
        ID_Cluster_mat <- as.data.frame(ID_Cluster_mat)
    
                                        #plot the edge
                                        #get the X-coordinate and y-coordinate of pies
        aa <- p$data

        desc <- y_union$Description[match(rownames(ID_Cluster_mat), y_union$ID)]
        i <- match(desc, aa$name) 

        ID_Cluster_mat$x <- aa$x[i]
        ID_Cluster_mat$y <- aa$y[i]
    
                                        #Change the radius value to fit the pie plot
        ID_Cluster_mat$radius <- sqrt(aa$size[i] / sum(aa$size))
                                        #ID_Cluster_mat$radius <- sqrt(aa$size / pi)
    
        x_loc1 <- min(ID_Cluster_mat$x)
        y_loc1 <- min(ID_Cluster_mat$y)
    ## x_loc2 <- min(ID_Cluster_mat$x)
    ## y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))
        if(nrow(ID_Cluster_mat) > 4) {    
        p + scatterpie::geom_scatterpie(aes(x=x,y=y,r=radius), data=ID_Cluster_mat,
                                    cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
        
            coord_equal()+
            geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
            scatterpie::geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
                                           labeller=function(x) round(sum(aa$size)*(x^2))) +
            labs(fill = "Cluster")
    ## annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")
	} else {
	    p + geom_node_point(aes_(color=~color, size=~size)) +
            geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
            scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
            scale_size(range=c(3, 8))
	}

    }
}

