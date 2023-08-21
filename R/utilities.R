##' automatically split barplot or dotplot into several facets
##' 
##'
##' @param by one of 'row' or 'column'
##' @param scales wether 'fixed' or 'free'
##' @param levels set facet levels
##' @return a ggplot object
##' @export
autofacet <- function(by = 'row', scales = "free", levels = NULL) {
    structure(list(by = by,
                scales = scales,
                levels = levels), 
            class = "autofacet")
}



# has_package <- function(pkg){
    # if (!requireNamespace(pkg, quietly  = TRUE)) {
        # stop(paste0(pkg, " needed for this function to work. Please install it."),
            # call. = FALSE)
    # }
# }



##' @method as.data.frame compareClusterResult
##' @export
as.data.frame.compareClusterResult <- function(x, ...) {
    as.data.frame(x@compareClusterResult, ...)
}


##' Prepare pie data for genes in cnetplot.
##' The function only works for compareClusterResult
##'
##' @importFrom DOSE geneID
##' @param y a data.frame converted from compareClusterResult
##' @return a data.frame
##' @noRd
prepare_pie_gene <- function(y) {
    gene_pie <- tibble::as_tibble(y[,c("Cluster", "Description", "geneID")])
    gene_pie$geneID <- strsplit(gene_pie$geneID, '/')
    gene_pie2 <- as.data.frame(tidyr::unnest(gene_pie, cols=geneID))
    gene_pie2 <- unique(gene_pie2)
    prepare_pie_data(gene_pie2, pie =  "equal", type = "gene")
}


##' Prepare pie data for categories in cnetplot/emapplot.
##' The function only works for compareClusterResult
##'
##' @param enrichDf a data.frame converted from compareClusterResult
##' @param pie proportion of clusters in the pie chart, one of 'equal' (default)
##' or 'Count'
##' @return a data.frame
##' @noRd
prepare_pie_category <- function(enrichDf, pie = "equal") {
    pie <- match.arg(pie, c("equal", "count", "Count"))
    if (pie == "count") pie <- "Count"

    pie_data <- enrichDf[,c("Cluster", "Description", "Count")]
    pie_data[,"Description"] <- as.character(pie_data[,"Description"])
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
##' @export
##' @examples
##' color_palette(c("red", "yellow", "green"))
##' @author guangchuang yu
color_palette <- function(colors) {
    # has_package("grDevices")
    grDevices::colorRampPalette(colors)(n = 299)
}


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

    if (x@readable && x@keytype != "SYMBOL") {
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

# fc_palette <- function(fc) {
    # if (all(fc > 0, na.rm=TRUE)) {
        # palette <- color_palette(c("blue", "red"))
    # } else if (all(fc < 0, na.rm=TRUE)) {
        # palette <- color_palette(c("green", "blue"))
    # } else {
        ## palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
    # }
    # return(palette)
# }

update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        if (inherits(x, 'list')) {
            showCategory <- showCategory[showCategory %in% names(x)]
        } else {
            showCategory <- intersect(showCategory, x$Description)
        }
        return(showCategory)
    }

    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    n <- showCategory
    if (inherits(x, 'list')) {
        nn <- length(x)
    } else {
        nn <- nrow(x)
    }
    if (nn < n) {
        n <- nn
    }

    return(n)
}

extract_geneSets <- function(x, n) {
    n <- update_n(x, n)

    if (inherits(x, 'list')) {
        geneSets <- x
    } else {
        geneSets <- geneInCategory(x) ## use core gene for gsea result
        y <- as.data.frame(x)
        geneSets <- geneSets[y$ID]
        names(geneSets) <- y$Description
    }
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
    # if (type == "dot") {
    #     if (by == "rowPercentage") {
    #         p <- ggplot(clProf.reshape.df,
    #                     aes_(x = x, y = ~Description, size = ~Percentage))
    #     } else if (by == "count") {
    #         p <- ggplot(clProf.reshape.df,
    #                     aes_(x = x, y = ~Description, size = ~Count))
    #     } else if (by == "geneRatio") {
    #         p <- ggplot(clProf.reshape.df,
    #                     aes_(x = x, y = ~Description, size = ~GeneRatio))
    #     } else {
    #         ## nothing here
    #     }
    #     p <- ggplot(clProf.reshape.df,
    #                 aes_(x = x, y = ~Description, size = by))
    #     if (any(colnames(clProf.reshape.df) == colorBy)) {
    #         p <- p +
    #             geom_point() +
    #             aes_string(color=colorBy) +
    #             scale_color_continuous(low="red", high="blue",
    #                                    guide=guide_colorbar(reverse=TRUE))
    #         ## scale_color_gradientn(guide=guide_colorbar(reverse=TRUE), colors = sig_palette)
    #     } else {
    #         p <- p + geom_point(colour="steelblue")
    #     }
    # }
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




##' Get the distance of the label
##'
##' @param dimension one of 1 and 2
##' @param label_location label_location
##' @noRd
get_label_diss <- function(dimension, label_location) {
    nn <- nrow(label_location)
    label_dis <- matrix(NA, nrow = nn, ncol = nn)
    colnames(label_dis) <- rownames(label_dis) <- label_location$label
    for (i in seq_len(nn - 1)) {
        for (j in (i + 1):nn) {
        label_dis[i ,j] <- label_location[i, dimension] -  label_location[j, dimension]
        }
    }
    label_diss <- reshape2::melt(label_dis)
    label_diss <- label_diss[label_diss[,1] != label_diss[,2], ]
    label_diss <- label_diss[!is.na(label_diss[,3]), ]
    label_diss[, 1] <- as.character(label_diss[, 1])
    label_diss[, 2] <- as.character(label_diss[, 2])
    return(label_diss)
}



# adjust_location <- function(label_location, x_adjust, y_adjust) {
    # label_diss_x <- get_label_diss(1, label_location)
    # label_diss_y <- get_label_diss(2, label_location)

    # label_diss_large <- which(abs(label_diss_y[, 3]) < y_adjust) %>%
        # intersect(which(label_diss_y[, 3] > 0)) %>%
        # intersect(which(abs(label_diss_x[, 3]) < x_adjust))

    # label_diss_small <- which(abs(label_diss_y[, 3]) < y_adjust) %>%
        # intersect(which(label_diss_y[, 3] < 0)) %>%
        # intersect(which(abs(label_diss_x[, 3]) < x_adjust))

    # label_location[label_diss_y[label_diss_large, 1], 2] <- label_location[label_diss_y[label_diss_large, 2], 2] + y_adjust
    # label_location[label_diss_y[label_diss_small, 1], 2] <- label_location[label_diss_y[label_diss_small, 2], 2] - y_adjust
    # return(label_location)
# }



#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `yulab.utils::str_wrap`
#' @noRd
#' @importFrom yulab.utils str_wrap
default_labeller <- function(n) {
    fun <- function(str){
        str <- gsub("_", " ", str)
        yulab.utils::str_wrap(str, n)
    }
    
    structure(fun, class = "labeller")
}

# from hadley wickham in "https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html"
#' Suppressing output
#'
#' @param x some code
#' @noRd
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


#' Get segment.size value for ggrepel
#' @param default default value of ggrepel.segment.size
#' @noRd
get_ggrepel_segsize <- function(default = 0.2) {
    getOption("ggrepel.segment.size", default = default)
}

#' Get warning message of changing parameter name
#' @param parameter old parameter name
#' @param params_df data frame with three columns: "original", "listname", and "present"
#' @noRd
get_param_change_message <- function(parameter, params_df) {
    paste0("Use '", params_df[parameter, "listname"], 
           " = list(", params_df[parameter, "present"], 
           " = your_value)' instead of '", params_df[parameter, "original"],
         "'.\n The ", params_df[parameter, "original"],
          " parameter will be removed in the next version.")
} 