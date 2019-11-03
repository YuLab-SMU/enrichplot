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

emap_graph_build <- function(n1,y1,geneSets1,layout1,color1) {
    if (n1 == 1) {
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- y1$Description
        V(g)$color <- "red"
        return(ggraph(g))
		} else {
        id <- y1[,"ID"]
        geneSets1 <- geneSets1[id]
        n1 <- nrow(y1) #
        w <- matrix(NA, nrow=n1, ncol=n1)
        colnames(w) <- rownames(w) <- y1$Description
        
        for (i in 1:n1) {
            for (j in i:n1) {
                w[i,j] <- overlap_ratio(geneSets1[id[i]], geneSets1[id[j]])
            }
        }
        
        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
        g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
        ## g <- delete.edges(g, E(g)[wd[,3] < 0.05])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y1$Description)))
        
        cnt <- sapply(geneSets1[idx], length)
        V(g)$size <- cnt
        
        colVar <- y1[idx, color1]
        V(g)$color <- colVar
    }


    p <- ggraph(g, layout=layout1)
    
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }
    return(p)
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
    p <- emap_graph_build(n1=n,y1=y,geneSets1=geneSets,layout1=layout,color1=color)
    if(n==1) {
        p + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name))
    } else {
        p + geom_node_point(aes_(color=~color, size=~size)) +
            geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
            scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
            ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
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
    p <- emap_graph_build(n1=n,y1=y_union,geneSets1=geneSets,layout1=layout,color1=color)


    ## then add the pie plot	
    ## Get the matrix data for the pie plot
    data_pie <- as.matrix(y[,1:2])
    ID_unique <- unique(data_pie[,2])
    Cluster_unique <- unique(data_pie[,1])
    
    
    ID_Cluster_mat <- matrix(0,length(ID_unique),length(Cluster_unique))
    rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    for(i in seq_len(nrow(data_pie))) {
        ID_Cluster_mat[data_pie[i,2],data_pie[i,1]] <- 1
    }
    
    ##use "region" as group of geom_scatterpie()
    ID_Cluster_mat <- as.data.frame(ID_Cluster_mat)
    ID_Cluster_mat$region <- factor(1:dim(ID_Cluster_mat)[1])
    
                                        #plot the edge
                                        #get the X-coordinate and y-coordinate of pies
    aa <- p$data

    desc <- y_union$Description[match(ID_unique, y_union$ID)]
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
    
    p + scatterpie::geom_scatterpie(aes(x=x,y=y, group=region,r=radius), data=ID_Cluster_mat,
                                    cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-4)],color=NA) +
        coord_equal()+
        scatterpie::geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
                                           labeller=function(x) round(sum(aa$size)*(x^2))) +
        labs(fill = "Cluster")
    ## annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")

}


