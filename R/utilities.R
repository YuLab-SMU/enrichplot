##' Prepare pie data for genes in cnetplot.
##' The function only works for compareClusterResult
##'
##' @param y a data.frame converted from compareClusterResult
##' @return a data.frame
##' @noRd 
prepare_pie_gene <- function(y) {
    gene_pie <- tibble::as_tibble(y[,c("Cluster", "ID", "geneID")])
    gene_pie$geneID <- strsplit(gene_pie$geneID, '/')
    gene_pie2 <- as.data.frame(tidyr::unnest(gene_pie, cols=geneID))
    gene_pie2 <- unique(gene_pie2)
    prepare_pie_data(gene_pie2, pie =  "equal", type = "gene")
}


##' Prepare pie data for categories in cnetplot/emapplot.
##' The function only works for compareClusterResult
##'
##' @param y a data.frame converted from compareClusterResult
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default) or 'Count'
##' @return a data.frame
##' @noRd
prepare_pie_category <- function(y, pie = "equal") {
    pie <- match.arg(pie, c("equal", "count", "Count"))
    if (pie == "count") pie <- "Count"
    
    pie_data <- y[,c("Cluster", "ID", "Count")]

    prepare_pie_data(pie_data, pie = pie)
}




prepare_pie_data <- function(pie_data, pie = "equal",type = "category") {
    if(type == "category"){
        ID_unique <- unique(pie_data[,2])
    } else {
        ID_unique <- unique(pie_data[,3])
    }
    
    Cluster_unique <- unique(pie_data[,1])
    ID_Cluster_mat <- matrix(0, nrow = length(ID_unique), ncol = length(Cluster_unique))
    rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
    if(pie == "Count") {
        for(i in seq_len(nrow(pie_data))) {
            ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- pie_data[i,3]
        }
        for(kk in seq_len(ncol(ID_Cluster_mat))) {
            ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
        }
        return(ID_Cluster_mat)
    }
    for(i in seq_len(nrow(pie_data))) {
        if(type == "category"){
            ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- 1
        } else {
            ID_Cluster_mat[pie_data[i,3],pie_data[i,1]] <- 1
    }
        
    }
    return(ID_Cluster_mat)
}


##' create color palette for continuous data
##'
##'
##' @title color_palette
##' @param colors colors of length >=2
##' @return color vector
##' @importFrom grDevices colorRampPalette
##' @export
##' @examples
##' color_palette(c("red", "yellow", "green"))
##' @author guangchuang yu
color_palette <- function(colors) colorRampPalette(colors)(n = 299)

sig_palette <- color_palette(c("red", "yellow", "blue"))

heatmap_palette <- color_palette(c("red", "yellow", "green"))

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

fc_readable <- function(x, foldChange = NULL) {
    if (is.null(foldChange))
        return(NULL)

    if(x@readable) {
        gid <- names(foldChange)
        if (is(x, 'gseaResult')) {
            ii <- gid %in% names(x@geneList)
        } else {
            ii <- gid %in% x@gene
        }
        gid[ii] <- x@gene2Symbol[gid[ii]]
        names(foldChange) <- gid
    }
    return(foldChange)
}

fc_palette <- function(fc) {
    if (all(fc > 0, na.rm=TRUE)) {
        palette <- color_palette(c("blue", "red"))
    } else if (all(fc < 0, na.rm=TRUE)) {
        palette <- color_palette(c("green", "blue"))
    } else {
        palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
    }
    return(palette)
}

update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        return(showCategory)
    }

    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    n <- showCategory
    if (nrow(x) < n) {
        n <- nrow(x)
    }

    return(n)
}

extract_geneSets <- function(x, n) {
    n <- update_n(x, n)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description
    if (is.numeric(n)) {
        return(geneSets[1:n])
    }
    return(geneSets[n]) ## if n is a vector of Description
}

##' Internal plot function for plotting compareClusterResult
##'
##'
##' @title plotting-clusterProfile
##' @param clProf.reshape.df data frame of compareCluster result
##' @param x x variable
##' @param type one of dot and bar
##' @param by one of percentage and count
##' @param title graph title
##' @param font.size graph font size
##' @param colorBy one of pvalue or p.adjust
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 scale_color_continuous
##' @importFrom ggplot2 guide_colorbar
##' @importFrom DOSE theme_dose
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
plotting.clusterProfile <- function(clProf.reshape.df,
                                    x = ~Cluster,
                                    type = "dot",
                                    colorBy = "p.adjust",
                                    by = "geneRatio",
                                    title="",
                                    font.size=12) {
    Description <- Percentage <- Count <- Cluster <- GeneRatio <- p.adjust <- pvalue <- NULL # to satisfy codetools
    if (type == "bar") {
        if (by == "percentage") {
            p <- ggplot(clProf.reshape.df,
                        aes(x=Description, y = Percentage, fill=Cluster))
        } else if (by == "count") {
            p <- ggplot(clProf.reshape.df,
                        aes(x=Description, y = Count, fill=Cluster))
        } else {

        }
        p <- p +
            geom_bar() +
                coord_flip()
    }
    if (type == "dot") {
        if (by == "rowPercentage") {
            p <- ggplot(clProf.reshape.df,
                        aes_(x = x, y = ~Description, size = ~Percentage))
        } else if (by == "count") {
            p <- ggplot(clProf.reshape.df,
                        aes_(x = x, y = ~Description, size = ~Count))
        } else if (by == "geneRatio") {
            p <- ggplot(clProf.reshape.df,
                        aes_(x = x, y = ~Description, size = ~GeneRatio))
        } else {
            ## nothing here
        }
        if (any(colnames(clProf.reshape.df) == colorBy)) {
            p <- p +
                geom_point() +
                aes_string(color=colorBy) +
                scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))
            ## scale_color_gradientn(guide=guide_colorbar(reverse=TRUE), colors = sig_palette)
        } else {
            p <- p + geom_point(colour="steelblue")
        }
    }
    p <- p + xlab("") + ylab("") + ggtitle(title) +
        theme_dose(font.size)
    ## theme(axis.text.x = element_text(colour="black", size=font.size, vjust = 1)) +
    ##     theme(axis.text.y = element_text(colour="black",
    ##           size=font.size, hjust = 1)) +
    ##               ggtitle(title)+theme_bw()
    ## p <- p + theme(axis.text.x = element_text(angle=angle.axis.x,
    ##                    hjust=hjust.axis.x,
    ##                    vjust=vjust.axis.x))
    return(p)
}
