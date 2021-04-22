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
##' @param reductFunction The function used for dimension reduction, 
##' e.g. 'stats::cmdscale' (the default), 'vegan::metaMDS', or 'vegan::metaMDS'.
##' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
##' @param ... additional parameters
##' 
##' additional parameters can refer the emapplot function: \link{emapplot}.
ssplot.enrichResult <- function(x, showCategory = 30, 
                                nCluster = NULL, reductFunction = "stats::cmdscale", 
                                with_edge = FALSE, ...) {      

    # arg <- match.call(expand.dots=TRUE)
    # arg <- as.list(environment())
    # reductFunctionName <- as.character(arg["reductFunction"])
    if (!is.character(reductFunction))
        stop("The value of reductFunction parameter must be a character")

    wrongMessage <- paste("Wrong reductFunction parameter or unsupported", 
         "dimensionality reduction method;",
         "set to default `reductFunction = 'stats::cmdscale'`")                     
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
    res <- reduction_dim(similarity_mat = termsim2, reductFunction = reductFunction, 
        wrongMessage = wrongMessage)
    coords <- get_coords(dim_reduction_data = res, wrongMessage = wrongMessage,
        similarity_mat  = termsim2, reductFunction = reductFunction)  
    p <- emapplot(x = x, showCategory = showCategory, 
                  coords = coords,
                  nCluster = nCluster,                 
                  with_edge = with_edge,
                  ...)

    ## Set axis label according to reductFunction
    p <- adj_axis(p = p, dim_reduction_data = res, reductFunction = reductFunction)
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
                                        reductFunction = "stats::cmdscale", 
                                        with_edge = FALSE, 
                                        cex_pie2axis = 0.0125, ...) { 
    # arg <- match.call(expand.dots=FALSE)
    # reductFunctionName <- as.character(arg["reductFunction"])
    wrongMessage <- paste("Wrong reductFunction parameter or unsupported", 
        "dimensionality reduction method;",
        "set to default `reductFunction = 'stats::cmdscale'`")   
    if (!is.character(reductFunction))
        stop("The value of reductFunction parameter must be a character")

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
    res <- reduction_dim(similarity_mat = termsim2, reductFunction = reductFunction, 
        wrongMessage = wrongMessage)
    coords <- get_coords(dim_reduction_data = res, wrongMessage = wrongMessage,
        similarity_mat  = termsim2, reductFunction = reductFunction)  
    p <- emapplot(x, showCategory = showCategory,
                  coords = coords,
                  split = split, pie = pie,
                  nCluster = nCluster,  
                  with_edge = with_edge, 
                  cex_pie2axis = cex_pie2axis, ...)
    ## Set axis label according to the method parameter
    p <- adj_axis(p = p, dim_reduction_data = res, reductFunction = reductFunction) 
    return(p + theme_classic())
}


##' Dimensionality reduction 
##'
##' @param similarity_mat Similarity matrix of terms.
##' @param reductFunction The function used for dimension reduction.
##' @param wrongMessage The information displayed when an error occurs.
##'
##' @noRd
reduction_dim <- function(similarity_mat, reductFunction, wrongMessage) {
    reductFunctionName <- gsub(".*::", "", reductFunction)   
    reductFunction <- eval(parse(text=reductFunction))
    ## If the similarity between the two terms is 1, 
    ## an error will be reported, so fine-tuning.
    similarity_mat[which(similarity_mat == 1)] <- 0.99999
    for (i in seq_len(nrow(similarity_mat))) similarity_mat[i, i] <- 1
    
    dune.dist <- stats::as.dist(1- similarity_mat)

    if (reductFunctionName == "cmdscale" || reductFunctionName == "wcmdscale") {
        return(reductFunction(dune.dist, eig = TRUE))   
    }

    if (reductFunctionName == "dudi.pco") {
        return(reductFunction(dune.dist, scannf = FALSE))
    }

    if (reductFunctionName == "sammon") {
        y <- stats::cmdscale(dune.dist)
        ## If the matrix y has duplicate rows it will report an error, so perturb slightly  
        dup <- which(duplicated(y) == TRUE)  
        y[dup, 1] <- y[dup, 1] + 10^-7 * seq_len(length(dup))  
        res <- reductFunction(d = dune.dist, y = y)
        return(res)
    }

    res <- tryCatch(expr  = reductFunction(dune.dist), 
                    error = function(e) {
                        message(wrongMessage)
                        stats::cmdscale(dune.dist, eig = TRUE)
                    })
    return(res)
}

# reduction_dim <- function(similarity_mat, nCluster, method) {
#     ## If the similarity between the two terms is 1, an error will be reported, so fine-tuning.
#     similarity_mat[which(similarity_mat == 1)] <- 0.99999
#     for (i in seq_len(nrow(similarity_mat))) similarity_mat[i, i] <- 1
    
#     dune.dist <- stats::as.dist(1- similarity_mat)
#     # res <- ape::pcoa(dune.dist, correction="none")
#     if (method == "cmdscale") {
#         res <- stats::cmdscale(dune.dist, eig = T)   
#     }

#     if (method == "sammon") {
#         y <- stats::cmdscale(dune.dist)
#         ## If the matrix y has duplicate rows it will report an error, so perturb slightly  
#         dup <- which(duplicated(y) == TRUE)  
#         y[dup, 1] <- y[dup, 1] + 10^-7 * seq_len(length(dup))  
#         res <- MASS::sammon(d = dune.dist, y = y)
#     }  
    
#     if (method == "monoMDS" || method == "isoMDS") {
#         res <- vegan::metaMDS(dune.dist, engine = method)
#     }
#     # coords <- as.data.frame(res$points[, 1:2])
#     # colnames(coords) <- c('x', 'y')

#     # if (method == "cmdscale") {
#     #     coords$pcoas <- 0
#     #     coords$pcoas[seq_len(length(res$eig))] <- as.numeric(res$eig)
#     # } else {
#     #     coords$stress <- res$stress
#     # }
#     return(res)
# }

##' Adjust axis label according to the dimension reduction method
##'
##' @param p ggplot2 object
##' @param dim_reduction_data a matrix data of dimension reduction result
##' @param reductFunction The name of function used for dimension reduction.
##' @noRd
adj_axis <- function(p, dim_reduction_data, reductFunction) {
    reductFunctionName <- gsub(".*::", "", reductFunction)   
    reductFunction <- eval(parse(text=reductFunction))
    pcoas <- NULL  
    eig <- dim_reduction_data$eig
    if (!is.null(eig) && length(eig) > 0) 
        pcoas <- as.numeric(eig)
    if (reductFunctionName == "pcoa") {
        pcoas <- as.numeric(dim_reduction_data$values$Eigenvalues)
    }
    title = NULL
    if (!is.null(pcoas)) {
        xlab = paste("Dimension1 (", format(100 * pcoas[1] / sum(pcoas), digits=4), "%)", sep = "")
        ylab = paste("Dimension2 (", format(100 * pcoas[2] / sum(pcoas), digits=4), "%)", sep = "")
    } else {
        xlab = "Dimension1"
        ylab = "Dimension2"
        if (!is.null(dim_reduction_data$stress)) {
            title = paste("stress = ", format(dim_reduction_data$stress, digits=4), sep = "")
        }       
    }
    p <- p +  labs(x = xlab, y = ylab, title = title) 
    return(p)
}

##' Extract the coordinates.
##'
##' Extract the coordinates of the first two dimensions from the dimensionality reduction result.
##'
##' @param dim_reduction_data a matrix data of dimension reduction result
##' @param wrongMessage The information displayed when an error occurs.
##' @param similarity_mat Similarity matrix of terms.
##' @param reductFunction The name of function used for dimension reduction.
##' @noRd
get_coords <- function(dim_reduction_data, 
                       wrongMessage, similarity_mat,
                       reductFunction) {
    reductFunctionName <- gsub(".*::", "", reductFunction)   
    reductFunction <- eval(parse(text=reductFunction))
    coords <- NULL
    if (reductFunctionName == "mds") {
        coords <- as.data.frame(dim_reduction_data$conf[, 1:2])
    } else if (reductFunctionName == "pcoa") {
        coords <- as.data.frame(dim_reduction_data$vectors[,1:2])
    } else if (reductFunctionName == "pco") {
        ## ecodist::pco
        if (!is.null(dim_reduction_data$vectors)) {
            coords <- as.data.frame(dim_reduction_data$vectors[, 1:2])
        }
        ## ade4::dudi.pco
        if (!is.null(dim_reduction_data$tab)) {
            coords <- as.data.frame(dim_reduction_data$tab[, 1:2])
        }
        ## labdsv::pco
        if (!is.null(dim_reduction_data$points)) {
            coords <- as.data.frame(dim_reduction_data$points[, 1:2])
        }
    } else if (reductFunctionName == "dudi.pco") {
        coords <- as.data.frame(dim_reduction_data$tab[, 1:2])
    } else {
        coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    }  
    if (is.null(coords) || nrow(coords) == 0) {
        message(wrongMessage)
        dim_reduction_data <- reduction_dim(similarity_mat = similarity_mat, reductFunction = stats::cmdscale, 
            wrongMessage = wrongMessage)
        coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    }
    colnames(coords) <- c('x', 'y')
    return(coords)
}
