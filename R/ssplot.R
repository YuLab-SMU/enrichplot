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
##' e.g. stats::cmdscale (the default), vegan::metaMDS, or ape::pcoa.
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
    # drfunName <- as.character(eval(substitute(alist(drfun))))
    wrongMessage <- paste("Wrong drfun parameter or unsupported",
         "dimensionality reduction method;",
         "set to default `drfun = 'stats::cmdscale'`")

    x <- quiet(dr(x = x, drfun = drfun, showCategory = showCategory))
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

    x <- quiet(dr(x = x, drfun = drfun,
                  showCategory = showCategory, 
                  split = split, pie = pie))
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



##' Dimension reduction
##'
##' @param x Enrichment result
##' @param drfun dimension reduction function 
##' @param showCategory A number or a vector of terms. If it is a number, 
##' the first n terms will be displayed. If it is a vector of terms, 
##' the selected terms will be displayed.
##' @param split separate result by 'category' variable
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) and 'Count'
dr <- function(x, drfun, showCategory, split = NULL, pie = NULL) {
    sim = get_pairwise_sim(x = x, showCategory = showCategory,
        split = split, pie = pie)
    ## If the similarity between the two terms is 1,
    ## an error will be reported in some method, so fine-tuning.
    sim[which(sim == 1)] <- 0.99999
    for (i in seq_len(nrow(sim))) sim[i, i] <- 1
    # drfun <- eval(parse(text=drfunName)) 
    # drfunName <- gsub(".*::", "", drfunName)
    drfun_env <- environment(fun = drfun) %>%
        rlang::env_name() %>%
        sub(".*:", "", .)
    distance_matrix <- stats::as.dist(1- sim)
    class(distance_matrix) <- c(drfun_env, class(distance_matrix)) 
    drs <- as.drs(distance_matrix = distance_matrix, drfun = drfun)
    x@dr <- drs
    return(x)
}


#' Dimensionality reduction and extract the results
#'
#' @param distance_matrix Distance matrix
#' @param drfun The function used for dimension reduction
#' @noRd
as.drs <- function(distance_matrix, drfun) UseMethod("as.drs")

as.drs.default <- function(distance_matrix, drfun) {
    wrongMessage <- paste("Wrong drfun parameter or unsupported",
        "dimensionality reduction method;",
        "set to default `drfun = 'stats::cmdscale'`")       
    message(wrongMessage)
    dim_reduction_data <- stats::cmdscale(distance_matrix, eig = TRUE)
    coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$eig)
    list(coords = coords, pcoa = pcoa)
}



#' @method as.drs stats
as.drs.stats <- function(distance_matrix, drfun) {
    dim_reduction_data <- drfun(distance_matrix, eig = TRUE)
    coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$eig)
    list(coords = coords, pcoa = pcoa)
}

#' @method as.drs MASS
as.drs.MASS <- function(distance_matrix, drfun) {
    y <- stats::cmdscale(distance_matrix)
    ## If the matrix y has duplicate rows it will report an error, so perturb slightly  
    dup <- which(duplicated(y) == TRUE)  
    y[dup, 1] <- y[dup, 1] + 10^-7 * seq_len(length(dup))  
    dim_reduction_data <- drfun(d = distance_matrix, y = y)
    coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    stress <- dim_reduction_data$stress
    list(coords = coords, stress = format(stress, digits=4))
}

#' @method as.drs vegan
as.drs.vegan <- function(distance_matrix, drfun) {
    if("engine" %in% names(as.list(rlang::get_expr(drfun)))) {
       dim_reduction_data <- drfun(distance_matrix, engine = "monoMDS")
       coords <- as.data.frame(dim_reduction_data$points[, 1:2])
       stress <- dim_reduction_data$stress
       return(list(coords = coords, stress = format(stress, digits=4)))
    } else {
       dim_reduction_data <- drfun(distance_matrix, eig = TRUE)
       coords <- as.data.frame(dim_reduction_data$points[, 1:2])
       pcoa <- as.numeric(dim_reduction_data$eig)
       return(list(coords = coords, pcoa = pcoa))
    }
}

#' @method as.drs ape
as.drs.ape <- function(distance_matrix, drfun) {
    dim_reduction_data <- drfun(distance_matrix)
    coords <- as.data.frame(dim_reduction_data$vectors[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$values$Eigenvalues)
    list(coords = coords, pcoa = pcoa)
}

#' @method as.drs smacof
as.drs.smacof<- function(distance_matrix, drfun) {
    dim_reduction_data <- drfun(distance_matrix)
    coords <- as.data.frame(dim_reduction_data$conf[, 1:2])
    stress <- dim_reduction_data$stress
    list(coords = coords, stress = format(stress, digits=4))
}


#' @method as.drs ecodist
as.drs.ecodist <- function(distance_matrix, drfun) {
    dim_reduction_data <- drfun(distance_matrix)
    coords <- as.data.frame(dim_reduction_data$vectors[, 1:2])
    rownames(coords) <- attr(distance_matrix, "Labels")
    pcoa <- as.numeric(dim_reduction_data$values)
    list(coords = coords, pcoa = pcoa)
}

#' @method as.drs labdsv
as.drs.labdsv <- function(distance_matrix, drfun) {
    dim_reduction_data <- drfun(distance_matrix)
    coords <- as.data.frame(dim_reduction_data$points[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$eig)
    list(coords = coords, pcoa = pcoa)
}



#' @method as.drs ade4
as.drs.ade4 <- function(distance_matrix, drfun) {
    dim_reduction_data <- drfun(distance_matrix, scannf = FALSE)
    coords <- as.data.frame(dim_reduction_data$tab[, 1:2])
    pcoa <- as.numeric(dim_reduction_data$eig)
    list(coords = coords, pcoa = pcoa)
}



get_pairwise_sim <- function(x, showCategory, split = NULL, pie = NULL) {
    if (class(x) == "compareClusterResult") {
        # y <- get_selected_category(showCategory, enrichResult, split)
        y <- fortify(model = x, showCategory = showCategory,
                     includeAll = TRUE, split = split)
        y$Cluster <- sub("\n.*", "", y$Cluster)
        keep <- rownames(prepare_pie_category(y, pie=pie))
    } else {
        n <- update_n(x, showCategory)
        if (is.numeric(n)) {
            keep <- seq_len(n)
        } else {
            keep <- match(n, rownames(x@termsim))
        }
    }
    if (length(keep) == 0) {
        stop("no enriched term found...")
    }
    fill_termsim(x, keep)
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

