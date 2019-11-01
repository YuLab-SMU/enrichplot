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
          function(x, showCategory = 5, color = "p.adjust", layout = "kk", by="geneRatio",split=NULL, includeAll=TRUE,...) {
              emapplot.compareClusterResult(x, showCategory = showCategory, color=color, layout = layout, by="geneRatio",split=split, includeAll=includeAll, ...)
          })
		  



graph_build <- function(n1,y1,geneSets1,layout1,color1) {
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
                w[i,j] = enrichplot:::overlap_ratio(geneSets1[id[i]], geneSets1[id[j]])
            }
        }

        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
        g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
		#g <- delete.edges(g, E(g)[wd[,3] < 0.05])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y1$Description)))

        cnt <- sapply(geneSets1[idx], length)
        V(g)$size <- cnt

        colVar <- y[idx, color1]
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
    p <- graph_build(n1=n,y1=y,geneSets1=geneSets,layout1=layout,color1=color)
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

# Make appropriate changes to several functions
geneID.compareClusterResult <- function(x) {
    as.character(x@compareClusterResult$geneID)
}
geneInCategory.compareClusterResult <- function(x) {
    setNames(strsplit(geneID(x), "/", fixed = TRUE), x@compareClusterResult$ID)

}
merge_compareClusterResult <- function(yy) {
	yy_union<- yy[!duplicated(yy$ID),]
	rownames(yy_union) <- yy_union$ID
	# Make a for loop to merge genes
	for(i in seq_len(dim(yy)[1])) {
	    gene1 <- unlist(strsplit(yy[i,"geneID"],"/", fixed = TRUE))
		gene2 <- unlist(strsplit(yy_union[yy[i,"ID"],"geneID"],"/", fixed = TRUE))
		genes <- union(gene1,gene2)
		gene_ID_uni <- paste(genes,collapse="/")
		yy_union[yy[i,"ID"],"geneID"] <- gene_ID_uni	
	}
	return(yy_union)
}


merge_compareClusterResult <- function(yy) {
    #one strsplit, one paste

	

	list_yy <- lapply(X = yy[,"geneID"],FUN=function(aa) {unlist(strsplit(as.character(aa),"/", fixed = TRUE))})
	yy_union<- yy[!duplicated(yy$ID),]
	list_yy_union <- lapply(X = yy_union[,"geneID"],FUN=function(aa) {unlist(strsplit(as.character(aa),"/", fixed = TRUE))})
	names(list_yy_union) <- yy_union$ID
	rownames(yy_union) <- yy_union$ID
	
	# Make a for loop to merge genes
	for(i in seq_len(dim(yy)[1])) {
	    #gene1 <- unlist(strsplit(yy[i,"geneID"],"/", fixed = TRUE))
		gene1 <- list_yy[[i]]
		#gene2 <- unlist(strsplit(yy_union[yy[i,"ID"],"geneID"],"/", fixed = TRUE))
		gene2 <- list_yy_union[yy[i,"ID"]][[1]]
		genes <- union(gene1,gene2)
		gene_ID_uni <- paste(genes,collapse="/")
		yy_union[yy[i,"ID"],"geneID"] <- gene_ID_uni	
	}
	return(yy_union)
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
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 coord_equal
##' @importFrom ggplot2 annotate
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @importFrom DOSE geneID
##' @method fortify compareClusterResult
##' @importFrom scatterpie geom_scatterpie
##' @importFrom stats setNames
emapplot.compareClusterResult <- function(x, showCategory = 5, color = "p.adjust", layout = "kk", by="geneRatio",split=NULL, includeAll=TRUE,...) {
region <- radius <- NULL

#pretreatment of x, just like dotplot do
y <- clusterProfiler:::fortify.compareClusterResult(x, showCategory=showCategory, by=by, includeAll=includeAll, split=split)
clusterr <- sapply(X = y$Cluster, FUN = function(dd){unlist(strsplit(as.character(dd),"\n"))[1]})
y$Cluster <- clusterr
geneSets <- geneInCategory(x) 

#geneSets <- geneInCategory(x) ## use core gene for gsea result

# Data structure transformation, combining the same ID (Description) genes
n <- update_n(x, showCategory)

y_union <- merge_compareClusterResult(y)

 if (n == 0) {
        stop("no enriched term found...")
    } 
	
p <- graph_build(n1=n,y1=y_union,geneSets1=geneSets,layout1=layout,color1=color)


# then add the pie plot	
#Get the matrix data for the pie plot
data_pie <- as.matrix(y[,1:2])
ID_unique <- unique(data_pie[,2])
Cluster_unique <- unique(data_pie[,1])


ID_Cluster_mat <- matrix(0,length(ID_unique),length(Cluster_unique))
rownames(ID_Cluster_mat) <- ID_unique
colnames(ID_Cluster_mat) <- Cluster_unique
for(i in seq_len(dim(data_pie)[1])) {
    ID_Cluster_mat[data_pie[i,2],data_pie[i,1]] <- ID_Cluster_mat[data_pie[i,2],data_pie[i,1]]+1
}

#use "region" as group of geom_scatterpie()
ID_Cluster_mat <- as.data.frame(ID_Cluster_mat)
ID_Cluster_mat$region <- factor(1:dim(ID_Cluster_mat)[1])

#plot the edge
#get the X-coordinate and y-coordinate of pies
aa <- p$data
	
ID_Cluster_mat$x <- aa$x
ID_Cluster_mat$y <- aa$y

#Change the radius value to fit the pie plot
ID_Cluster_mat$radius <- sqrt(aa$size / sum(aa$size))
#ID_Cluster_mat$radius <- sqrt(aa$size / pi)

x_loc1 <- min(ID_Cluster_mat$x)
y_loc1 <- min(ID_Cluster_mat$y)
x_loc2 <- min(ID_Cluster_mat$x)
y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))

p + geom_scatterpie(aes(x=x,y=y, group=region,r=radius), data=ID_Cluster_mat,
                           cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-4)],color=NA) + coord_equal()+
						    scatterpie::geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,labeller=function(x) round(sum(aa$size)*(x^2)))+
							 annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")



}


