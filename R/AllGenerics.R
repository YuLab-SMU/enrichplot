##' Gene-Concept Network
##'
##'
##' plot linkages of genes and enriched concepts
##' (e.g. GO categories, KEGG pathways)
##' @title cnetplot
##' @rdname cnetplot
##' @param x enrichment result
##' @param showCategory number of enriched terms to display
##' @param foldChange fold Change
##' @param layout layout of the network
##' @param ... additional parameters
##' @return ggplot object
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' cnetplot(x)
cnetplot <- function(x, showCategory = 5,
                     foldChange = NULL, layout = "kk", ...) {
    UseMethod("cnetplot", x)
}

##' dotplot for enrichment result
##'
##'
##' @title dotplot
##' @rdname dotplot
##' @param height input object, just name it to make it consistent with barplot
##' @param x variable for x-axis, one of 'geneRatio' or 'Count'
##' @param color variable that used to color enriched terms,
##'              e.g. pvalue, p.adjust or qvalue
##' @param showCategory number of enriched terms to display
##' @param split separate result by 'category' variable
##' @param font.size font size
##' @param title plot title
##' @return plot
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' dotplot(x)
dotplot <- function(height, x = "geneRatio", color = "p.adjust",
                    showCategory=10, split = NULL, font.size=12, title = "") {
    UseMethod("dotplot", height)
}

##' Enrichment Map for enrichment result of
##' over-representation test or gene set enrichment analysis
##'
##'
##' This function visualizes gene sets as a network (i.e. enrichment map).
##' Mutually overlapping gene sets tend to cluster together, making it easier for interpretation.
##' @title emapplot
##' @rdname emapplot
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue, p.adjust or qvalue
##' @param layout layout of the map
##' @param ... additional parameter
##' @return ggplot object
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' emapplot(x)
emapplot <- function(x, showCategory = 30, color="p.adjust", layout = "kk", ...) {
    UseMethod("emapplot", x)
}

##' plot induced GO DAG of significant terms
##'
##'
##' @title goplot
##' @rdname goplot
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue, p.adjust or qvalue
##' @param layout layout of the map
##' @param geom label geom, one of 'label' or 'text'
##' @param ... additional parameter
##' @return ggplot object
##' @export
goplot <- function(x, showCategory = 10, color = "p.adjust", layout = "sugiyama", geom = "text", ...) {
    UseMethod("goplot", x)
}

##' visualize analyzing result of GSEA
##'
##' plotting function for gseaResult
##' @title gseaplot
##' @rdname gseaplot
##' @param x object of gsea result
##' @param geneSetID geneSet ID
##' @param by one of "runningScore" or "position"
##' @param title plot title
##' @param ... additional parameters
##' @return ggplot2 object
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' x <- gseDO(geneList)
##' gseaplot(x, geneSetID=1)
gseaplot <- function(x, geneSetID, by = "all", title = "", ...) {
    UseMethod("gseaplot", x)
}

##' heatmap like plot for functional classification
##'
##'
##' @title heatplot
##' @rdname heatplot
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param foldChange fold Change
##' @export
##' @return ggplot object
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' heatplot(x)
##' @author guangchuang yu
heatplot <- function(x, showCategory=30, foldChange=NULL) {
    UseMethod("heatplot", x)
}


##' ridgeline plot for GSEA result
##'
##'
##' @title ridgeplot
##' @rdname ridgeplot
##' @param x gseaResult object
##' @param showCategory number of categories for plotting
##' @param fill one of "pvalue", "p.adjust", "qvalue"
##' @param core_enrichment whether only using core_enriched genes
##' @return ggplot object
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' x <- gseDO(geneList)
##' ridgeplot(x)
ridgeplot <- function(x, showCategory=30, fill="p.adjust", core_enrichment = TRUE) {
    UseMethod("ridgeplot", x)
}


#' geneInCategory generic
#'
#' @param x enrichResult
#' @return 'geneInCategory' return a list of genes, by spliting the input gene vector to enriched functional categories
#' @export
#' @examples
#' library(DOSE)
#' data(geneList)
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' geneInCategory(x)
geneInCategory <- function(x) {
   UseMethod("geneInCategory", x)
}


#' geneID generic
#'
#' @param x enrichResult object
#' @return 'geneID' return the 'geneID' column of the enriched result which can be converted to data.frame via 'as.data.frame'
#' @export
#' @examples
#' library(DOSE)
#' data(geneList)
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' geneID(x)
geneID <- function(x) {
   UseMethod("geneID", x)
}
