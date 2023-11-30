##' Gene-Concept Network
##'
##'
##' plot linkages of genes and enriched concepts
##' (e.g. GO categories, KEGG pathways)
##' @title cnetplot
##' @rdname cnetplot
##' @param x Enrichment result.
##' @param showCategory A number or a vector of terms. If it is a number, 
##' the first n terms will be displayed. If it is a vector of terms, 
##' the selected terms will be displayed.
##' @param foldChange Fold Change of nodes, the default value is NULL. 
##' If the user provides the Fold Change value of the nodes, 
##' it can be used to set the color of the gene node.
##' Will be removed in the next version.
##' @param layout Layout of the map, e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
##' 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
##' @param ... Additional parameters
##' @return ggplot object
##' @export
##' @examples
##' \dontrun{
##'     library(DOSE)
##'     data(geneList)
##'     de <- names(geneList)[1:100]
##'     x <- enrichDO(de)
##'     x2 <- pairwise_termsim(x)
##'     cnetplot(x2)
##'     # use `layout` to change the layout of map
##'     cnetplot(x2, layout = "star")
##'     # use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
##'     cnetplot(x2, showCategory = 10)
##'     categorys <- c("pre-malignant neoplasm", "intestinal disease",
##'                    "breast ductal carcinoma", "non-small cell lung carcinoma")
##'     cnetplot(x2, showCategory = categorys)
##'     # 'compareClusterResult' object is also supported.
##'     library(clusterProfiler)
##'     library(DOSE)
##'     library(org.Hs.eg.db)
##'     data(gcSample)
##'     xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
##'     xx2 <- pairwise_termsim(xx)
##'     cnetplot(xx2)
##' }
setGeneric("cnetplot",
           function(x,  ...)
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
##' \dontrun{
##'     library(DOSE)
##'     data(geneList)
##'     de <- names(geneList)[1:100]
##'     x <- enrichDO(de)
##'     dotplot(x)
##'     # use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
##'     dotplot(x, showCategory = 10)
##'     categorys <- c("pre-malignant neoplasm", "intestinal disease",
##'                    "breast ductal carcinoma", "non-small cell lung carcinoma")
##'     dotplot(x, showCategory = categorys)
##'     # It can also graph compareClusterResult
##'     data(gcSample)
##'     library(clusterProfiler)
##'     library(DOSE)
##'     library(org.Hs.eg.db)
##'     data(gcSample)
##'     xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
##'     xx2 <- pairwise_termsim(xx)
##'     library(ggstar)
##'     dotplot(xx2)
##'     dotplot(xx2, shape = TRUE)
##'     dotplot(xx2, group = TRUE)
##'     dotplot(xx2, x = "GeneRatio", group = TRUE, size = "count")
##' }
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
##' easier for interpretation. When the similarity between terms meets 
##' a certain threshold (default is 0.2, adjusted by parameter `min_edge`),
##' there will be edges between terms. The stronger the similarity, 
##' the shorter and thicker the edges. The similarity between terms is 
##' obtained by function `pairwise_termsim`, the details of similarity 
##' calculation can be found in its documentation: \link{pairwise_termsim}.
##' @title emapplot
##' @rdname emapplot
##' @param x Enrichment result.
##' @param showCategory A number or a vector of terms. If it is a number, 
##' the first n terms will be displayed. If it is a vector of terms, 
##' the selected terms will be displayed.
##' @param ... Additional parameters
##' @return ggplot object
##' @export
##' @examples
##' \dontrun{
##'     library(DOSE)
##'     data(geneList)
##'     de <- names(geneList)[1:100]
##'     x <- enrichDO(de)
##'     x2 <- pairwise_termsim(x)
##'     emapplot(x2)
##'     # use `layout` to change the layout of map
##'     emapplot(x2, layout = "star")
##'     # use `showCategory` to  select the displayed terms. It can be a number of a vector of terms.
##'     emapplot(x2, showCategory = 10)
##'     categorys <- c("pre-malignant neoplasm", "intestinal disease",
##'                    "breast ductal carcinoma")
##'     emapplot(x2, showCategory = categorys)
##' 
##'     # It can also graph compareClusterResult
##'     library(clusterProfiler)
##'     library(DOSE)
##'     library(org.Hs.eg.db)
##'     data(gcSample)
##'     xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
##'     xx2 <- pairwise_termsim(xx)
##'     emapplot(xx2)
##' }
setGeneric("emapplot",
           function(x,  ...)
               standardGeneric("emapplot")
           )



           
##' Get the similarity matrix
##'
##'
##' This function add similarity matrix to the termsim slot of enrichment result.
##' Users can use the `method` parameter to select the method of calculating similarity.
##' The Jaccard correlation coefficient(JC) is used by default, and it applies to all situations.
##' When users want to calculate the correlation between GO terms or DO terms, they can also choose
##' "Resnik", "Lin", "Rel" or "Jiang" (they are semantic similarity calculation methods from GOSemSim packages),
##' and at this time, the user needs to provide `semData` parameter, which can be obtained through
##' \link{godata} function in GOSemSim package.
##' @title pairwise_termsim
##' @rdname pairwise_termsim
##' @param x enrichment result.
##' @param method method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC"(Jaccard similarity coefficient) methods.
##' @param semData GOSemSimDATA object, can be obtained through
##' \link{godata} function in GOSemSim package.
##' @param showCategory number of enriched terms to display, default value is 200.
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
           function(x, method = "JC", semData = NULL, showCategory = 200)
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
##' @examples
##' \dontrun{
##' 	library(clusterProfiler)
##'   data(geneList, package = "DOSE")
##' 	de <- names(geneList)[1:100]
##' 	yy <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
##'     goplot(yy)
##'     goplot(yy, showCategory = 5)
##' }
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
##' @param foldChange fold Change.
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param ... Additional parameters
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
           function(x, showCategory = 30, ...)
               standardGeneric("heatplot")
           )



##' ridgeline plot for GSEA result
##'
##'
##' @title ridgeplot
##' @rdname ridgeplot
##' @param x gseaResult object
##' @param showCategory A number or a vector of terms. If it is a number, 
##' the first n terms will be displayed. If it is a vector of terms, 
##' the selected terms will be displayed.
##' @param fill one of "pvalue", "p.adjust", "qvalue"
##' @param core_enrichment whether only using core_enriched genes
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param ... additional parameters
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
                    label_format = 30, ...)
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


##' Functional grouping tree diagram for enrichment result of 
##' over-representation test or gene set enrichment analysis.
##'
##'
##' This function visualizes gene sets as a tree.
##' Gene sets with high similarity tend to cluster together, making it easier
##' for interpretation.
##' @title treeplot
##' @rdname treeplot
##' @param x enrichment result.
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue,
##' p.adjust or qvalue
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param ... additional parameters
##' @return ggplot object
##' @export
##' @examples
##' \dontrun{
##'     library(clusterProfiler)
##'     library(org.Hs.eg.db)
##'     library(enrichplot)
##'     library(GOSemSim)
##'     library(ggplot2)
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
##'     ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
##'     treeplot(ego2, showCategory = 30)
##'     # use `hilight = FALSE` to remove ggtree::geom_hilight() layer.
##'     treeplot(ego2, showCategory = 30, hilight = FALSE)
##'     # use `offset` parameter to adjust the distance of bar and tree.
##'     treeplot(ego2, showCategory = 30, hilight = FALSE, offset = rel(1.5))
##'     # use `offset_tiplab` parameter to adjust the distance of nodes and branches.
##'     treeplot(ego2, showCategory = 30, hilight = FALSE, offset_tiplab = rel(1.5))
##'     keep <- rownames(ego2@termsim)[c(1:10, 16:20)]
##'     keep
##'     treeplot(ego2, showCategory = keep)
##'     treeplot(ego2, showCategory = 20, 
##'         group_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"))
##'     # It can also graph compareClusterResult
##'     data(gcSample)
##'     xx <- compareCluster(gcSample, fun="enrichKEGG",
##'                          organism="hsa", pvalueCutoff=0.05)
##'     xx <- pairwise_termsim(xx)                     
##'     treeplot(xx)                     
##'     
##'     # use `geneClusterPanel` to change the gene cluster panel.
##'     treeplot(xx, geneClusterPanel = "dotplot")  
##'     
##'     treeplot(xx, geneClusterPanel = "pie")  
##'    }
setGeneric("treeplot",
           function(x, ...)
               standardGeneric("treeplot")
           )
           
##' Similarity space plot of enrichment analysis results.
##'
##' @title ssplot
##' @rdname ssplot
##' @inheritParams emapplot
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
##'         ont           = "BP",
##'         pAdjustMethod = "BH",
##'         pvalueCutoff  = 0.01,
##'         qvalueCutoff  = 0.05,
##'         readable      = TRUE)
##'     d <- godata('org.Hs.eg.db', ont="BP")
##'     ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
##'     ssplot(ego2)    
##' }
setGeneric("ssplot",
           function(x, ...)
               standardGeneric("ssplot")
           )
           
                      
