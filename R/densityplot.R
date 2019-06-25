##' plot logFC distribution of selected gene sets
##'
##' 
##' @title gseadist 
##' @param x GSEA result
##' @param IDs gene set IDs
##' @param type one of 'density' or 'boxplot'
##' @return distribution plot
##' @importFrom ggplot2 geom_density
##' @importFrom ggplot2 geom_boxplot
##' @export
##' @author Guangchuang Yu
gseadist <- function(x, IDs, type =  'density') {
    d <- data.frame(gene = names(x@geneList),
                    logFC = x@geneList,
                    category = 'All Genes')

    ds <- do.call('rbind', lapply(IDs, function(i) {
        if (!is.numeric(i)) {
            i <- match(i, x$ID)
            if (is.na(i))
                i <- match(i,  x$Description)
        } 
        id <- x$ID[i]

        gene <- x@geneSets[[id]]
        gs <- x@geneList[gene]
        gs <- gs[!is.na(gs)]
        data.frame(gene = names(gs),
                   logFC = gs,
                   category = x$Description[i])
    }))
    dd <- rbind(d, ds)

    p <- ggplot(dd) + theme_minimal()

    if (type == 'density') {
        p <- p + 
            geom_density(aes_(x = ~logFC, color = ~category)) +
            ## geom_rug(data = ds, show.legend = FALSE) +
            ylab(NULL) +
            theme(legend.title = element_blank(),
                  legend.position = 'bottom')
    } else if (type == 'boxplot') {
        p <- p +
            geom_boxplot(aes_(x = ~category, y = ~logFC, fill = ~category))  +
            xlab(NULL) +
            theme(legend.position = 'none')
    } 
    return(p)
}
