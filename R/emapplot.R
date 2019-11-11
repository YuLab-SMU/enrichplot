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




##' Get graph.data.frame() result
##'
##' @param y a data.frame of clusterProfiler result
##' @param geneSets a list gene sets with the names of enrichment IDs
##' @param color a string, the column name of y for nodes colours
##' @return result of graph.data.frame()
##' @noRd
emap_graph_build <- function(y,geneSets,color) {
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

        for (i in seq_len(n-1)) {
            for (j in (i+1):n) {
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
##' @param pie_scale scale of pie plot 
##' @author Guangchuang Yu
emapplot.enrichResult <- function(x, showCategory = 30, color="p.adjust", layout = "kk", pie_scale = 1,...) {
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

    g <- emap_graph_build(y=y,geneSets=geneSets,color=color)
    if(n == 1) {
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    }
    ##} else {

    p <- ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }
    p + geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8) * pie_scale)

}



##' Merge the compareClusterResult file
##'
##' @param yy a data.frame of clusterProfiler result
##'
##' @return a data.frame
##' @noRd
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

##' Prepare the data for the pie plot
##'
##' @param y a data.frame of clusterProfiler result
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @return a data.frame
##' @noRd
deal_data_pie <- function(y, pie = "equal") {
    pie <- match.arg(pie, c("equal", "count", "Count"))
    if (pie == "count") pie <- "Count"

    data_pie <- as.matrix(y[,c(1,2,10)])
    ##rownames(data_pie) <- paste0(data_pie[,1],data_pie[,2])
    ID_unique <- unique(data_pie[,2])
    Cluster_unique <- unique(data_pie[,1])
    ID_Cluster_mat <- matrix(0,length(ID_unique),length(Cluster_unique))
    rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
    if(pie == "Count") {
        for(i in seq_len(nrow(data_pie))) {
            ID_Cluster_mat[data_pie[i,2],data_pie[i,1]] <- data_pie[i,3]
        }
        for(kk in seq_len(ncol(ID_Cluster_mat))) {
            ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
        }
        return(ID_Cluster_mat)
    }
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
##' @importFrom scatterpie geom_scatterpie
##' @importFrom scatterpie geom_scatterpie_legend
##' @importFrom DOSE geneInCategory
##' @importFrom DOSE geneID
##' @importClassesFrom DOSE compareClusterResult
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @param legend_n number of circle in legend
##' @param pie_scale scale of pie plot 
##' @method fortify compareClusterResult
##' @importFrom scatterpie geom_scatterpie
##' @importFrom stats setNames
##' @noRd
emapplot.compareClusterResult <- function(x, showCategory = 5, color = "p.adjust",
                                          layout = "kk", split=NULL, pie = "equal",
                                          legend_n = 5, pie_scale = 1, ...) {

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
    g <- emap_graph_build(y=y_union,geneSets=geneSets,color=color)
    ## when y just have one line
    if(is.null(dim(y)) | nrow(y) == 1) {
        title <- y$Cluster
        p <- ggraph(g) + geom_node_point(color="red", size=5 * pie_scale) +
            geom_node_text(aes_(label=~name)) + theme_void() +
            labs(title=title)
        return(p)
    }

    if(is.null(dim(y_union)) | nrow(y_union) == 1) {
         ##return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
        p <- ggraph(g)
        ID_Cluster_mat <- deal_data_pie(y, pie=pie)

        ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*pie_scale)
        colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],"x","y","radius")


        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                                 cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA)+
            xlim(-3,3) + ylim(-3,3) + coord_equal()+ geom_node_text(aes_(label=~name), repel=TRUE) +
            theme_void()+labs(fill = "Cluster")
        return(p)

    }
    p <- ggraph(g, layout=layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }

    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- deal_data_pie(y,pie=pie)


                                        #plot the edge
                                        #get the X-coordinate and y-coordinate of pies
    aa <- p$data

    desc <- y_union$Description[match(rownames(ID_Cluster_mat), y_union$ID)]
    i <- match(desc, aa$name)

    ID_Cluster_mat$x <- aa$x[i]
    ID_Cluster_mat$y <- aa$y[i]

                                        #Change the radius value to fit the pie plot
    ID_Cluster_mat$radius <- sqrt(aa$size[i] / sum(aa$size)) * pie_scale
                                        #ID_Cluster_mat$radius <- sqrt(aa$size / pi)

    x_loc1 <- min(ID_Cluster_mat$x)
    y_loc1 <- min(ID_Cluster_mat$y)
    ## x_loc2 <- min(ID_Cluster_mat$x)
    ## y_loc2 <- min(ID_Cluster_mat$y)+0.1*(max(ID_Cluster_mat$y)-min(ID_Cluster_mat$y))
    if(ncol(ID_Cluster_mat) > 4) {
        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                                 cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +

            coord_equal()+
            geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
            geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1, n = legend_n,
                                   labeller=function(x) round(sum(aa$size)*((x/pie_scale)^2))) +
            labs(fill = "Cluster")
        return(p)
    }
    ## annotate("text", label = "gene number", x = x_loc2, y = y_loc2, size = 4, colour = "red")
    title <- colnames(ID_Cluster_mat)[1]
    p + geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8) * pie_scale)  +labs(title= title)
}

