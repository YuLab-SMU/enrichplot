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
##' @param method The function used for MDS, 
##' one of "cmdscale" (the default), "sammon", "monoMDS" and "isoMDS".
##' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
##' @param ... additional parameters
##' 
##' additional parameters can refer the emapplot function: \link{emapplot}.
ssplot.enrichResult <- function(x, showCategory = 30, 
                                nCluster = NULL, method = "cmdscale", 
                                with_edge = FALSE, ...) {
    n <- update_n(x, showCategory)
    ## Use the similarity matrix to reduce the dimension and redistribute the node coordinates.
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

    ## Dimensionality reduction and K-means clustering using similarity matrix
    data <- get_cluster(termsim2, nCluster = nCluster, method = method)

    p <- emapplot(x = x, showCategory = showCategory, 
                  layout = "ssplot",
                  nCluster = nCluster,
                  data = data,
                  with_edge = with_edge,
                  ...)

    ## Set axis label according to method
    p <- adj_axis(p = p, data = data)
    return(p + theme_classic())
}




##' @rdname ssplot
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 coord_equal
##' @importClassesFrom DOSE compareClusterResult
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'
##' @importFrom stats setNames
ssplot.compareClusterResult <- function(x, showCategory = 30,
                                        split = NULL, pie = "equal",
                                        nCluster = NULL, 
                                        method = "cmdscale", 
                                        with_edge = FALSE, ...) {

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
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    termsim2 <- fill_termsim(x, rownames(ID_Cluster_mat))

    ## Dimensionality reduction and K-means clustering using similarity matrix
    data <- get_cluster(termsim2, nCluster = nCluster, method = method)

    p <- emapplot(x, showCategory = 30,
                  layout = "nicely",
                  split = NULL, pie = "equal",
                  nCluster = NULL, data = data, 
                  with_edge = with_edge, ...)
    ## Set axis label according to the method parameter
    p <- adj_axis(p = p, data = data)
   
    return(p + theme_classic())
}


##' Group terms
##'
##' @param termsim2 Similarity matrix of terms.
##' @param nCluster Numeric, the number of cluster.
##' @param method The function used for MDS, 
##' one of "cmdscale" (the default), "sammon", "monoMDS" and "isoMDS".
##'
##' @noRd
get_cluster <- function(termsim2, nCluster, method) {
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

    if (method == "cmdscale") {
        data$pocas <- 0
        data$pocas[seq_len(length(res$eig))] <- as.numeric(res$eig)
    } else {
        data$stress <- res$stress
    }
    return(data)
}
