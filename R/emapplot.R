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
          function(x, showCategory = 30, color = "p.adjust", layout = "kk", ...) {
              emapplot.compareClusterResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })
		  
		  
		  
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
    } else if (n == 1) {
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- y$Description
        V(g)$color <- "red"
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    } else {
        id <- y[,1]
        geneSets <- geneSets[id]

        n <- nrow(y) #
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$Description

        for (i in 1:n) {
            for (j in i:n) {
                w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }

        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
        g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))

        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt

        colVar <- y[idx, color]
        V(g)$color <- colVar
    }


    p <- ggraph(g, layout=layout)

    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }

    p + geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8))
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

geom_scatterpie_legend2 <- function (radius, x, y, n = 5, labeller) 
{
    if (length(radius) > n) {
        radius <- unique(sapply(seq(min(radius), max(radius), 
            length.out = n), scatterpie:::round_digit))
    }
    label <- FALSE
    if (!missing(labeller)) {
        if (!inherits(labeller, "function")) {
            stop("labeller should be a function for converting radius")
        }
        label <- TRUE
    }
    dd <- data.frame(r = radius, start = 0, end = 2 * pi, x = x, 
        y = y + radius - max(radius), maxr = max(radius))
    if (label) {
        dd$label <- labeller(dd$r)
    }
    else {
	# Change the label to the actual number of genes
        dd$label <- (dd$r)^2*500
    }
    list(ggforce::geom_arc_bar(aes_(x0 = ~x, y0 = ~y, r0 = ~r, r = ~r, 
        start = ~start, end = ~end), data = dd, inherit.aes = FALSE), 
        geom_segment(aes_(x = ~x, xend = ~x + maxr * 1.5, y = ~y + 
            r, yend = ~y + r), data = dd, inherit.aes = FALSE), 
        geom_text(aes_(x = ~x + maxr * 1.6, y = ~y + r, label = ~label), 
            data = dd, hjust = "left", inherit.aes = FALSE))
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
##' @importFrom scatterpie geom_scatterpie
##' @importFrom stats setNames
emapplot.compareClusterResult <- function(x, showCategory = 30, color="p.adjust", layout = "kk", ...) {
region <- radius <- NULL
n <- enrichplot:::update_n(x, showCategory)
geneSets <- geneInCategory(x) 
y <- as.data.frame(x)
  if (is.numeric(n)) {
        y <- y[1:n,]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }

#geneSets <- geneInCategory(x) ## use core gene for gsea result

# Data structure transformation, combining the same ID (Description) genes


y_union <- merge_compareClusterResult(y)

    if (n <= 1) {
        stop("no enriched term found...")
    } else {
	    id <- y_union$ID
		geneSets <- geneSets[id]
        n <- nrow(y_union) #
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y_union$Description

        for (i in 1:n) {
            for (j in i:n) {
                w[i,j] = enrichplot:::overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }

        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
		g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.05])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y_union$Description)))

        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt

        colVar <- y_union[idx, color]
        V(g)$color <- colVar
	}
	

# then add the pie plot	
data_pie <- as.matrix(y[,1:2])
ID_unique <- unique(data_pie[,2])
Cluster_unique <- unique(data_pie[,1])
ID_Cluster_mat <- matrix(0,length(ID_unique),length(Cluster_unique))
rownames(ID_Cluster_mat) <- ID_unique
colnames(ID_Cluster_mat) <- Cluster_unique
for(i in seq_len(dim(data_pie)[1])) {
    ID_Cluster_mat[data_pie[i,2],data_pie[i,1]] <- ID_Cluster_mat[data_pie[i,2],data_pie[i,1]]+1
}
ID_Cluster_mat <- as.data.frame(ID_Cluster_mat)
ID_Cluster_mat$region <- factor(1:dim(ID_Cluster_mat)[1])

  p <- ggraph(g, layout=layout)

    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }

aa <- p$data
	
ID_Cluster_mat$x <- aa$x
ID_Cluster_mat$y <- aa$y


ID_Cluster_mat$radius <- sqrt(aa$size / 500)

x_loc1 <- min(ID_Cluster_mat$x)
y_loc1 <- min(ID_Cluster_mat$y)
x_loc2 <- min(ID_Cluster_mat$x)
y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))

p + geom_scatterpie(aes(x=x,y=y, group=region,r=radius), data=ID_Cluster_mat,
                           cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-4)],color=NA) + coord_equal()+
						    geom_scatterpie_legend2(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1)+
							 annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")
}


