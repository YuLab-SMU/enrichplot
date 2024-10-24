#' cnetplot
#' 
#' category-gene-network plot
#' @rdname cnetplot
#' @param x input object
#' @param layout network layout
#' @param showCategory selected category to be displayed
#' @param color_category color of category node
#' @param size_category relative size of the category
#' @param color_item color of item node
#' @param size_item relative size of the item (e.g., genes)
#' @param color_edge color of edge
#' @param size_edge relative size of edge
#' @param node_label one of 'all', 'none', 'category' and 'item'
#' @param foldChange numeric values to color the item (e.g, foldChange of gene expression values)
#' @param hilight selected category to be highlighted
#' @param hilight_alpha transparent value for not selected to be highlight
#' @param ... additional parameters
#' @importFrom ggtangle cnetplot
#' @method cnetplot enrichResult
#' @export
#' @seealso
#' [cnetplot][ggtangle::cnetplot]
cnetplot.enrichResult <- function(
        x, layout = igraph::layout_with_kk,
        showCategory = 5,
        color_category= "#E5C494", size_category = 1, 
        color_item = "#B3B3B3", size_item = 1, 
        color_edge = "grey", size_edge=.5,
        node_label = "all", 
        foldChange = NULL,
        hilight = "none",
        hilight_alpha = .3,
        ...) {

    geneSets <- extract_geneSets(x, showCategory)
    foldChange <- fc_readable(x, foldChange)    

    p <- cnetplot(geneSets, layout = layout, showCategory = showCategory, foldChange = foldChange, ...)
    p <- p + set_enrichplot_color(colors = get_enrichplot_color(3), name = "fold change")
    if (!is.null(foldChange)) {
        p <- p + guides(size  = guide_legend(order = 1), 
                        color = guide_colorbar(order = 2))
    }

    return(p + guides(alpha = "none"))
}

#' @rdname cnetplot
#' @method cnetplot gseaResult
#' @export
cnetplot.gseaResult <- cnetplot.enrichResult

#' @rdname cnetplot
#' @param pie one of 'equal' or 'Count' to set the slice ratio of the pies
#' @method cnetplot compareClusterResult
#' @export
cnetplot.compareClusterResult <- function(
        x, layout = igraph::layout_with_kk,
        showCategory = 5,
        color_category= "#E5C494", size_category = 1, 
        color_item = "#B3B3B3", size_item = 1, 
        color_edge = "grey", size_edge=.5,
        node_label = "all", 
        foldChange = NULL,
        hilight = "none",
        hilight_alpha = .3,
        pie = "equal",
        ...) {

    d <- tidy_compareCluster(x, showCategory)
    y <- split(d$geneID, d$Description)
    gs <- lapply(y, function(item) unique(unlist(strsplit(item, split="/"))))

    p <- cnetplot(gs, layout = layout, showCategory=length(gs), foldChange = foldChange, size_category=0, ...)

    add_node_pie(p, d, pie)
}

#' @importFrom ggplot2 coord_fixed
add_node_pie <- function(p, d, pie = "equal") {
    dd <- d[,c('Cluster', 'Description', 'Count')]
    if (pie == "equal") dd$Count <- 1
    dd <- tidyr::pivot_wider(dd, names_from=.data$Cluster, values_from=.data$Count, values_fill=0)
    
    p <- p %<+% dd +
        scatterpie::geom_scatterpie(cols=as.character(unique(d$Cluster)), legend_name = "Cluster", color=NA) +
        coord_fixed() +
        guides(size = "none")

    return(p)
}

tidy_compareCluster <- function(x, showCategory) {
    d <- fortify(x, showCategory = showCategory, includeAll = TRUE, split = NULL)
    d$Cluster <- sub("\n.*", "", d$Cluster)

    if ("core_enrichment" %in% colnames(d)) { ## for GSEA result
        d$geneID <- d$core_enrichment
    }
    return(d)
}
