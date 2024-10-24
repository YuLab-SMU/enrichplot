#' @method cnetplot enrichResult
#' @export
#' @seealso
#' [cnetplot][ggtangle::cnetplot]
cnetplot.enrichResult <- function(x, showCategory = 5, foldChange = NULL, ...) {
    geneSets <- extract_geneSets(x, showCategory)
    foldChange <- fc_readable(x, foldChange)    

    p <- cnetplot(geneSets, showCategory = showCategory, foldChange = foldChange, ...)
    p <- p + set_enrichplot_color(colors = get_enrichplot_color(3), name = "fold change")
    if (!is.null(foldChange)) {
        p <- p + guides(size  = guide_legend(order = 1), 
                        color = guide_colorbar(order = 2))
    }

    return(p + guides(alpha = "none"))
}

#' @method cnetplot gseaResult
#' @export
cnetplot.gseaResult <- cnetplot.enrichResult

#' @param pie one of 'equal' or 'Count' to set the slice ratio of the pies
#' @method cnetplot compareClusterResult
#' @export
cnetplot.compareClusterResult <- function(x, showCategory = 5, foldChange = NULL, pie = "equal", ...) {
    d <- x@compareClusterResult

    y <- split(d$geneID, d$Description)
    gs <- lapply(y, function(item) unique(unlist(strsplit(item, split="/"))))

    p <- cnetplot(gs, showCategory=showCategory, foldChange = foldChange, ...)

    dd <- d[,c('Cluster', 'Description', 'Count')]

    if (pie == "equal") dd$Count <- 1

    dd <- tidyr::pivot_wider(dd, names_from=Cluster, values_from=Count, values_fill=0)

    p <- p %<+% yy +
        scatterpie::geom_scatterpie(cols=as.character(unique(d$Cluster)), legend_name = "Cluster") +
        coord_fixed() +
        guides(size = "none")
    return(p)
}

