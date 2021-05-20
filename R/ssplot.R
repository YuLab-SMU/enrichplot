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
##' @param drfun The function used for dimension reduction,
##' e.g. 'stats::cmdscale' (the default), 'vegan::metaMDS', or 'vegan::metaMDS'.
##' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
##' @param ... additional parameters
##'
##' additional parameters can refer the emapplot function: \link{emapplot}.
ssplot.enrichResult <- function(x, showCategory = 30,
                                nCluster = NULL, drfun = stats::cmdscale,
                                with_edge = FALSE, ...) {
    if (is.character(drfun)) {
        drfun <- eval(parse(text=drfun))
    }  
    drfunName <- as.character(eval(substitute(alist(drfun))))
    wrongMessage <- paste("Wrong drfun parameter or unsupported",
         "dimensionality reduction method;",
         "set to default `drfun = 'stats::cmdscale'`")

    x <- dr(enrichResult = x, drfun = drfun, drfunName = drfunName, showCategory = showCategory)
    coords <- x@dr$coords
    colnames(coords) <- c("x", "y")
    p <- emapplot(x = x, showCategory = showCategory,
                  coords = coords,
                  nCluster = nCluster,
                  with_edge = with_edge,
                  ...)

    ## Set axis label according to drfun
    p <- adj_axis(p = p, drs = x@dr)
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
                                        drfun = stats::cmdscale,
                                        with_edge = FALSE,
                                        cex_pie2axis = 0.0125, ...) {
    if (is.character(drfun)) {
        drfun <- eval(parse(text=drfun))
    }  
    drfunName <- as.character(eval(substitute(alist(drfun))))
    ## Use the similarity matrix to reduce the dimension and redistribute the node coordinates.
    ## Dimensionality reduction
    x <- dr(enrichResult = x, drfun = drfun, drfunName = drfunName, 
        showCategory = showCategory, split = split, pie = pie)
    coords <- x@dr$coords
    colnames(coords) <- c("x", "y")
    p <- emapplot(x, showCategory = showCategory,
                  coords = coords,
                  split = split, pie = pie,
                  nCluster = nCluster,
                  with_edge = with_edge,
                  cex_pie2axis = cex_pie2axis, ...)
    ## Set axis label according to the method parameter
    p <- adj_axis(p = p, drs = x@dr)
    return(p + theme_classic())
}



dr <- function(enrichResult, drfun, drfunName, showCategory, split = NULL, pie = NULL) {

    sim = get_pairwise_sim(enrichResult = enrichResult, showCategory = showCategory,
        split = split, pie = pie)   
    ## If the similarity between the two terms is 1,
    ## an error will be reported in some method, so fine-tuning.
    sim[which(sim == 1)] <- 0.99999
    for (i in seq_len(nrow(sim))) sim[i, i] <- 1

    drfunName <- gsub(".*::", "", drfunName)
    distance_matrix <- stats::as.dist(1- sim)
    class(distance_matrix) <- c(drfunName, class(distance_matrix))
    result <- reduction_dim(distance_matrix = distance_matrix, drfun = drfun)
    class(result) <- c(drfunName, class(result))
    drs <- as.drs(result)
    # drs$method <- drfunName
    enrichResult@dr <- drs
    return(enrichResult)
}


#' as.drs
#'
#' @param distance_matrix Distance matrix
#' @param ... extra args
#' @noRd
reduction_dim <- function(distance_matrix, ...) UseMethod("reduction_dim")

reduction_dim.default <- function(distance_matrix, drfun) {
    wrongMessage <- paste("Wrong drfun parameter or unsupported",
        "dimensionality reduction method;",
        "set to default `drfun = 'stats::cmdscale'`")
    res <- tryCatch(expr  = drfun(distance_matrix,),
                error = function(e) {
                    message(wrongMessage)
                    stats::cmdscale(distance_matrix, eig = TRUE)
                })
}

#' @method reduction_dim cmdscale
reduction_dim.cmdscale <- function(distance_matrix, drfun) {
    drfun(distance_matrix, eig = TRUE)
}

#' @method reduction_dim wcmdscale
reduction_dim.wcmdscale <- function(distance_matrix, drfun) {
    drfun(distance_matrix, eig = TRUE)
}

#' @method reduction_dim dudi.pco
reduction_dim.dudi.pco <- function(distance_matrix, drfun) {
    drfun(distance_matrix, scannf = FALSE)
}

#' @method reduction_dim sammon
reduction_dim.sammon <- function(distance_matrix, drfun) {
    y <- stats::cmdscale(distance_matrix)
    ## If the matrix y has duplicate rows it will report an error, so perturb slightly
    dup <- which(duplicated(y) == TRUE)
    y[dup, 1] <- y[dup, 1] + 10^-7 * seq_len(length(dup))
    drfun(d = distance_matrix, y = y)
}


get_pairwise_sim <- function(enrichResult, showCategory, split = NULL, pie = NULL) {
    if (class(enrichResult) == "compareClusterResult") {
        y <- get_selected_category(showCategory, enrichResult, split)
        keep <- rownames(prepare_pie_category(y, pie=pie))
    } else {
        n <- update_n(enrichResult, showCategory)
        if (is.numeric(n)) {
            keep <- seq_len(n)
        } else {
            keep <- match(n, rownames(enrichResult@termsim))
        }   
    }  
    if (length(keep) == 0) {
        stop("no enriched term found...")
    }  
    fill_termsim(enrichResult, keep)
}


#' as.drs
#'
#' @param dim_reduction_data dim reduction data
#' @noRd
as.drs <- function(dim_reduction_data) UseMethod("as.drs")

as.drs.default <- function(dim_reduction_data) {
    coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$eig)
    list(coords = coords, pcoa = pcoa)
}

#' @method as.drs mds
as.drs.mds <- function(dim_reduction_data) {
    coords <- as.data.frame(dim_reduction_data$conf[, 1:2])
    stress <- dim_reduction_data$stress
    list(coords = coords, stress = format(dim_reduction_data$stress, digits=4))
}

#' @method as.drs pcoa
as.drs.pcoa <- function(dim_reduction_data) {
    coords <- as.data.frame(dim_reduction_data$vectors[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$values$Eigenvalues)
    list(coords = coords, pcoa = pcoa)
}

#' @method as.drs pco
as.drs.pco <- function(dim_reduction_data) {
    if (!is.null(dim_reduction_data$vectors)) {
        coords <- as.data.frame(dim_reduction_data$vectors[, 1:2])
    }
    ## labdsv::pco
    if (!is.null(dim_reduction_data$points)) {
        coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    }
    pcoa <- as.numeric(dim_reduction_data$eig)
    if (is.null(pcoa) || length(pcoa) == 0) {
        pcoa <- as.numeric(dim_reduction_data$values)
    }
    list(coords = coords, pcoa = pcoa)
}

#' @method as.drs dudi.pco
as.drs.dudi.pco <- function(dim_reduction_data) {
    coords <- as.data.frame(dim_reduction_data$tab[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$eig)
    list(coords = coords, pcoa = pcoa)
}

#' @method as.drs sammon
as.drs.sammon <- function(dim_reduction_data) {
    coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    stress <- dim_reduction_data$stress
    list(coords = coords, stress = format(dim_reduction_data$stress, digits=4))
}

#' @method as.drs metaMDS
as.drs.metaMDS <- function(dim_reduction_data) {
    as.drs.sammon(dim_reduction_data)
}


##' Adjust axis label according to the dimension reduction method
##'
##' @param p ggplot2 object
##' @param drs dimension reduction result
##' @noRd
adj_axis <- function(p, drs) {
    title = NULL
    pcoas <- drs$pcoa
    if (!is.null(pcoas)) {
        xlab = paste("Dimension1 (", format(100 * pcoas[1] / sum(pcoas), digits=4), "%)", sep = "")
        ylab = paste("Dimension2 (", format(100 * pcoas[2] / sum(pcoas), digits=4), "%)", sep = "")
    } else {
        xlab = "Dimension1"
        ylab = "Dimension2"
        if (!is.null(drs$stress)) {
            title = paste("stress = ", format(drs$stress, digits=4), sep = "")
        }
    }
    p <- p +  labs(x = xlab, y = ylab, title = title)
    return(p)
}






