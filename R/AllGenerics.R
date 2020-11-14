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
setGeneric("cnetplot",
           function(x, showCategory = 5,
                    foldChange = NULL, layout = "kk", ...)
               standardGeneric("cnetplot")
           )


##' dotplot for enrichment result
##'
##'
##' @title dotplot
##' @rdname dotplot
##' @param object input object
##' @param ... additional parameters
##' @return plot
##' @importFrom methods setGeneric
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' dotplot(x)
setGeneric("dotplot",
           function(object,  ...)
               standardGeneric("dotplot")
           )

##' Enrichment Map for enrichment result of
##' over-representation test or gene set enrichment analysis
##'
##'
##' This function visualizes gene sets as a network (i.e. enrichment map).
##' Mutually overlapping gene sets tend to cluster together, making it
##' easier for interpretation.
##' @title emapplot
##' @rdname emapplot
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue,
##' p.adjust or qvalue
##' @param layout layout of the map
##' @param ... additional parameters
##' @return ggplot object
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' de <- names(geneList)[1:100]
##' x <- enrichDO(de)
##' x2 <- pairwise_termsim(x)
##' emapplot(x2)
setGeneric("emapplot",
           function(x, showCategory = 30, color="p.adjust", layout = "kk", ...)
               standardGeneric("emapplot")
           )



##' Functional grouping network diagram for enrichment result of
##' over-representation test or gene set enrichment analysis
##'
##'
##' This function visualizes gene sets as a grouped network (i.e. enrichment map).
##' Gene sets with high similarity tend to cluster together, making it easier
##' for interpretation.
##' @title emapplot_cluster
##' @rdname emapplot_cluster
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue,
##' p.adjust or qvalue
##' @param ... additional parameters
##' @return ggplot object
##' @export
##' @examples
##' \dontrun{
##'     library(clusterProfiler)
##'     library(org.Hs.eg.db)
##'     library(enrichplot)
##'     library(GOSemSim)
##'     library(DOSE)
##'     data(geneList)
##'     gene <- names(geneList)[abs(geneList) > 2]
##'     ego <- enrichGO(gene  = gene,
##'         universe      = names(geneList),
##'         OrgDb         = org.Hs.eg.db,
##'         ont           = "CC",
##'         pAdjustMethod = "BH",
##'         pvalueCutoff  = 0.01,
##'         qvalueCutoff  = 0.05,
##'         readable      = TRUE)
##'     d <- godata('org.Hs.eg.db', ont="BP")
##'     ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
##'     emapplot_cluster(ego2, showCategory = 80)
##'    }
setGeneric("emapplot_cluster",
           function(x, showCategory = nrow(x), color="p.adjust", label_format = 30, ...)
               standardGeneric("emapplot_cluster")
           )

##' Get the similarity matrix
##'
##'
##' This function add similarity matrix to the termsim slot of enrichment result.
##' @title pairwise_termsim
##' @rdname pairwise_termsim
##' @param x enrichment result.
##' @param method method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC"(Jaccard similarity coefficient) methods.
##' @param semData GOSemSimDATA object
##' @param showCategory number of enriched terms to display
##' @examples
##' \dontrun{
##'     library(clusterProfiler)
##'     library(org.Hs.eg.db)
##'     library(enrichplot)
##'     library(GOSemSim)
##'     library(DOSE)
##'     data(geneList)
##'     gene <- names(geneList)[abs(geneList) > 2]
##'     ego <- enrichGO(gene  = gene,
##'         universe      = names(geneList),
##'         OrgDb         = org.Hs.eg.db,
##'         ont           = "BP",
##'         pAdjustMethod = "BH",
##'         pvalueCutoff  = 0.01,
##'         qvalueCutoff  = 0.05,
##'         readable      = TRUE)
##'     d <- godata('org.Hs.eg.db', ont="BP")
##'     ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
##'     emapplot(ego2)
##'     emapplot_cluster(ego2)
##'    }
setGeneric("pairwise_termsim",
           function(x, method = "JC", semData = NULL, showCategory = 30)
               standardGeneric("pairwise_termsim")
           )

##' plot induced GO DAG of significant terms
##'
##'
##' @title goplot
##' @rdname goplot
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue,
##' p.adjust or qvalue
##' @param layout layout of the map
##' @param geom label geom, one of 'label' or 'text'
##' @param ... additional parameter
##' @return ggplot object
##' @export
setGeneric("goplot",
           function(x, showCategory = 10, color = "p.adjust",
                    layout = "sugiyama", geom = "text", ...)
               standardGeneric("goplot")
           )

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
setGeneric("gseaplot",
           function(x, geneSetID, by = "all", title = "", ...) {
               standardGeneric("gseaplot")
           })


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
setGeneric("heatplot",
           function(x, showCategory=30, foldChange=NULL)
               standardGeneric("heatplot")
           )



##' ridgeline plot for GSEA result
##'
##'
##' @title ridgeplot
##' @rdname ridgeplot
##' @param x gseaResult object
##' @param showCategory number of categories for plotting
##' @param fill one of "pvalue", "p.adjust", "qvalue"
##' @param core_enrichment whether only using core_enriched genes
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' by default wraps names longer that 30 characters
##' @return ggplot object
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' x <- gseDO(geneList)
##' ridgeplot(x)
setGeneric("ridgeplot",
           function(x, showCategory=30, fill="p.adjust", core_enrichment = TRUE,
                    label_format = 30)
               standardGeneric("ridgeplot")
           )


##' upsetplot method generics
##'
##'
##' @docType methods
##' @name upsetplot
##' @rdname upsetplot-methods
##' @title upsetplot method
##' @param x object
##' @param ... additional parameters
##' @return plot
##' @export
setGeneric("upsetplot", function(x, ...) standardGeneric("upsetplot"))

