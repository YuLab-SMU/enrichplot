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
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' @param method The function used for dimension reduction, 
##' one of "cmdscale" (the default), "sammon", "monoMDS" and "isoMDS".
##' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
##' @param ... additional parameters
##' 
##' additional parameters can refer the emapplot function: \link{emapplot}.
ssplot.enrichResult <- function(x, showCategory = 30, 
                                nCluster = NULL, method = "cmdscale", 
                                with_edge = FALSE, ...) {
    method <- match.arg(method, c("cmdscale", "sammon", "monoMDS", "isoMDS"))                                
    n <- update_n(x, showCategory)
    if (is.numeric(n)) {
        keep <- seq_len(n)
    } else {
        keep <- match(n, rownames(x@termsim))
    }
    if (length(keep) == 0) {
        stop("no enriched term found...")
    }
    ## Complete similarity matrix
    termsim2 <- fill_termsim(x, keep)
    
    ## Dimensionality reduction 
    res <- reduction_dim(termsim2, nCluster = nCluster, method = method)
    coords <- as.data.frame(res$points[, 1:2])
    colnames(coords) <- c('x', 'y')
    p <- emapplot(x = x, showCategory = showCategory, 
                  coords = coords,
                  nCluster = nCluster,                 
                  with_edge = with_edge,
                  ...)

    ## Set axis label according to method
    p <- adj_axis(p = p, dim_reduction_data = res, method = method)
    return(p + theme_classic())
}




##' @rdname ssplot
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 coord_equal
##' @importClassesFrom DOSE compareClusterResult
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'
##' @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis, 
##' the default value is 0.0125.
##' @importFrom stats setNames
ssplot.compareClusterResult <- function(x, showCategory = 30,
                                        split = NULL, pie = "equal",
                                        nCluster = NULL, 
                                        method = "cmdscale", 
                                        with_edge = FALSE, 
                                        cex_pie2axis = 0.0125, ...) {
    method <- match.arg(method, c("cmdscale", "sammon", "monoMDS", "isoMDS"))  
    ## Complete the similarity matrix and retain the required nodes                                   
    y <- get_selected_category(showCategory, x, split)
    ## Then the similarity matrix is used to reduce the 
    ## dimension and redistribute the node coordinates.
    ## Complete similarity matrix
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    termsim2 <- fill_termsim(x, rownames(ID_Cluster_mat))

    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- merge_compareClusterResult(y)
     

    ## Use the similarity matrix to reduce the dimension and redistribute the node coordinates.
    ## Dimensionality reduction 
    res <- reduction_dim(termsim2, nCluster = nCluster, method = method)
    coords <- as.data.frame(res$points[, 1:2])
    colnames(coords) <- c('x', 'y')
    p <- emapplot(x, showCategory = showCategory,
                  coords = coords,
                  split = split, pie = pie,
                  nCluster = nCluster,  
                  with_edge = with_edge, 
                  cex_pie2axis = cex_pie2axis, ...)
    ## Set axis label according to the method parameter
    p <- adj_axis(p = p, dim_reduction_data = res, method = method)
   
    return(p + theme_classic())
}


##' Dimensionality reduction 
##'
##' @param similarity_mat Similarity matrix of terms.
##' @param nCluster Numeric, the number of cluster.
##' @param method The function used for dimension reduction, 
##' one of "cmdscale" (the default), "sammon", "monoMDS" and "isoMDS".
##'
##' @noRd
reduction_dim <- function(similarity_mat, nCluster, method) {
    ## If the similarity between the two terms is 1, an error will be reported, so fine-tuning.
    similarity_mat[which(similarity_mat == 1)] <- 0.99999
    for (i in seq_len(nrow(similarity_mat))) similarity_mat[i, i] <- 1
    
    dune.dist <- stats::as.dist(1- similarity_mat)
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
        res <- vegan::metaMDS(dune.dist, engine = method)
    }
    # coords <- as.data.frame(res$points[, 1:2])
    # colnames(coords) <- c('x', 'y')

    # if (method == "cmdscale") {
    #     coords$pocas <- 0
    #     coords$pocas[seq_len(length(res$eig))] <- as.numeric(res$eig)
    # } else {
    #     coords$stress <- res$stress
    # }
    return(res)
}

##' Adjust axis label according to the dimension reduction method
##'
##' @param p ggplot2 object
##' @param dim_reduction_data a matrix data of dimension reduction result
##' @param method The function used for dimension reduction, 
##' one of "cmdscale" (the default), "sammon", "monoMDS" and "isoMDS".
##' @noRd
adj_axis <- function(p, dim_reduction_data, method) {
    if (method == "cmdscale") {
        pocas <- as.numeric(dim_reduction_data$eig)
        xlab = paste("pcoa 1 (", format(100 * pocas[1] / sum(pocas), digits=4), "%)", sep = "")
        ylab = paste("pcoa 2 (", format(100 * pocas[2] / sum(pocas), digits=4), "%)", sep = "")
        title = "PCoA"
    } else {
        xlab = "Dimension1"
        ylab = "Dimension2"
        title = paste(method, " (stress = ", format(dim_reduction_data$stress, digits=4), ")", sep = "")
    }
    p <- p +  labs(x = xlab, y = ylab, title = title) 
    return(p)
}
