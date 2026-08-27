#' Dot plot for enrichment result
#'
#'
#' @title dotplot
#' @rdname dotplot
#' @param object input object.
#' @param ... additional parameters.
#' @return plot.
#' @importFrom methods setGeneric
#' @export
#' @examples
#' \dontrun{
#'     library(DOSE)
#'     data(geneList)
#'     de <- names(geneList)[1:100]
#'     x <- enrichDO(de)
#'     dotplot(x)
#'     # use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
#'     dotplot(x, showCategory = 10)
#'     categories <- c("pre-malignant neoplasm", "intestinal disease",
#'                    "breast ductal carcinoma", "non-small cell lung carcinoma")
#'     dotplot(x, showCategory = categories)
#'     # It can also graph compareClusterResult
#'     data(gcSample)
#'     library(clusterProfiler)
#'     library(DOSE)
#'     library(org.Hs.eg.db)
#'     data(gcSample)
#'     xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
#'     xx2 <- pairwise_termsim(xx)
#'     library(ggstar)
#'     dotplot(xx2)
#'     dotplot(xx2, shape = TRUE)
#'     dotplot(xx2, group = TRUE)
#'     dotplot(xx2, x = "GeneRatio", group = TRUE, size = "count")
#' }
#' @author Guangchuang Yu
setGeneric("dotplot", function(object, ...) {
    standardGeneric("dotplot")
})

#' Shared term-plot parameters
#'
#' @name enrichplot-term-params
#' @param showCategory number of categories to display, or a vector of terms.
#' @param color variable used to color enriched terms, e.g. `'pvalue'`,
#'   `'p.adjust'`, or `'qvalue'`.
#' @param label_format a numeric wrap width, or a custom function to format
#'   axis labels.
#' @keywords internal
NULL

#' Shared parameters for enrichment plots
#'
#' @name enrichplot-common-params
#' @param color variable used to color enriched terms, e.g. `'pvalue'`,
#'   `'p.adjust'`, or `'qvalue'`.
#' @param showCategory number of categories to display, or a vector of terms.
#' @param size variable used to scale category size, one of `"geneRatio"`,
#'   `"Percentage"`, or `"count"`.
#' @param split apply `showCategory` to each category specified by `split`,
#'   e.g., `"ONTOLOGY"`, `"category"`, or `"intersect"`. Default is `NULL`.
#' @param font.size font size.
#' @param title figure title.
#' @param label_format a numeric wrap width, or a custom function to format
#'   axis labels.
#' @param includeAll logical value.
#' @keywords internal
NULL

#' Enrichment Map for enrichment result of
#' over-representation test or gene set enrichment analysis
#'
#'
#' This function visualizes gene sets as a network (i.e. enrichment map).
#' Mutually overlapping gene sets tend to cluster together, making it
#' easier for interpretation. When the similarity between terms meets
#' a certain threshold (default is 0.2, adjusted by parameter `min_edge`),
#' there will be edges between terms. The stronger the similarity,
#' the shorter and thicker the edges. The similarity between terms is
#' obtained by the function `pairwise_termsim`. Details of the similarity
#' calculation can be found in its documentation: [pairwise_termsim()].
#' @title emapplot
#' @rdname emapplot
#' @param x Enrichment result.
#' @param showCategory number of categories to display or a vector of terms.
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'     library(DOSE)
#'     data(geneList)
#'     de <- names(geneList)[1:100]
#'     x <- enrichDO(de)
#'     x2 <- pairwise_termsim(x)
#'     emapplot(x2)
#'     # use `layout` to change the layout of map
#'     emapplot(x2, layout = "star")
#'     # use `showCategory` to  select the displayed terms. It can be a number of a vector of terms.
#'     emapplot(x2, showCategory = 10)
#'     categories <- c("pre-malignant neoplasm", "intestinal disease",
#'                    "breast ductal carcinoma")
#'     emapplot(x2, showCategory = categories)
#'
#'     # It can also graph compareClusterResult
#'     library(clusterProfiler)
#'     library(DOSE)
#'     library(org.Hs.eg.db)
#'     data(gcSample)
#'     xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
#'     xx2 <- pairwise_termsim(xx)
#'     emapplot(xx2)
#' }
#' @author Guangchuang Yu
setGeneric(
    "emapplot",
    function(
        x,
        layout = igraph::layout_with_kk,
        coords = NULL,
        showCategory = 30,
        color = "p.adjust",
        size_category = 1,
        min_edge = .2,
        color_edge = "grey",
        size_edge = .5,
        node_label = "category",
        node_label_size = 5,
        pie = "equal",
        layer = NULL,
        label_format = 30,
        clusterFunction = stats::kmeans,
        nWords = 4,
        nCluster = NULL,
        ...
    ) {
        standardGeneric("emapplot")
    }
)


#' Get the similarity matrix
#'
#'
#' This function adds a similarity matrix to the termsim slot of the enrichment result.
#' Users can use the `method` parameter to select the method of calculating the similarity.
#' The Jaccard correlation coefficient (JC) is used by default, and it applies to all situations.
#' When users want to calculate the correlation between GO terms or DO terms, they can also choose
#' "Resnik", "Lin", "Rel" or "Jiang" (they are semantic similarity calculation methods from the 'GOSemSim' package),
#' and at this time, the user needs to provide the `semData` parameter, which can be obtained through
#' [GOSemSim::godata()].
#' @title pairwise_termsim
#' @rdname pairwise_termsim
#' @param x enrichment result.
#' @param method method of calculating the similarity between nodes,
#' one of "Resnik", "Lin", "Rel", "Jiang", "Wang", and
#' "JC" (Jaccard similarity coefficient) methods.
#' @param semData `GOSemSimDATA` object, can be obtained through
#' `GOSemSim::godata`.
#' @param showCategory number of enriched terms to be calculated. The default value is the number of enriched terms, or 200 if the number of enriched terms exceeds 200.
#' @examples
#' \dontrun{
#'     library(clusterProfiler)
#'     library(org.Hs.eg.db)
#'     library(enrichplot)
#'     library(GOSemSim)
#'     library(DOSE)
#'     data(geneList)
#'     gene <- names(geneList)[abs(geneList) > 2]
#'     ego <- enrichGO(gene  = gene,
#'         universe      = names(geneList),
#'         OrgDb         = org.Hs.eg.db,
#'         ont           = "BP",
#'         pAdjustMethod = "BH",
#'         pvalueCutoff  = 0.01,
#'         qvalueCutoff  = 0.05,
#'         readable      = TRUE)
#'     d <- godata('org.Hs.eg.db', ont="BP")
#'     ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
#'     emapplot(ego2)
#'     emapplot_cluster(ego2)
#'    }
setGeneric(
    "pairwise_termsim",
    function(x, method = "JC", semData = NULL, showCategory = NULL) {
        standardGeneric("pairwise_termsim")
    }
)

#' Plot induced GO DAG of significant terms
#'
#'
#' @title goplot
#' @rdname goplot
#' @param x enrichment result.
#' @inheritParams enrichplot-term-params
#' @param layout layout of the map
#' @param geom label geom, one of 'label' or 'text'
#' @param ... additional parameters.
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' 	library(clusterProfiler)
#'   data(geneList, package = "DOSE")
#' 	de <- names(geneList)[1:100]
#' 	yy <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
#'     goplot(yy)
#'     goplot(yy, showCategory = 5)
#' }
#' @author Guangchuang Yu
setGeneric(
    "goplot",
    function(
        x,
        showCategory = 10,
        color = "p.adjust",
        layout = "sugiyama",
        geom = "text",
        ...
    ) {
        standardGeneric("goplot")
    }
)

#' Visualize GSEA analysis results
#'
#' Plotting function for gseaResult
#' @title gseaplot
#' @rdname gseaplot
#' @param x gseaResult object
#' @param geneSetID geneSet ID
#' @param by one of "runningScore" or "position"
#' @param title plot title
#' @param ... additional parameters
#' @return ggplot2 object
#' @export
#' @examples
#' \donttest{
#' library(DOSE)
#' data(geneList)
#' x <- gseDO(geneList)
#' gseaplot(x, geneSetID=1)
#' }
#' @author Guangchuang Yu
setGeneric("gseaplot", function(x, geneSetID, by = "all", title = "", ...) {
    standardGeneric("gseaplot")
})


#' Heatmap-like plot for functional classification
#'
#'
#' @title heatplot
#' @rdname heatplot
#' @param x enrichment result.
#' @param showCategory number of enriched terms to display
#' @param foldChange fold change.
#' @param label_format a numeric value setting the wrap length, alternatively a
#' custom function to format axis labels.
#' @param ... Additional parameters
#' @export
#' @return ggplot object
#' @examples
#' \dontrun{
#' library(DOSE)
#' data(geneList)
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' heatplot(x)
#' }
#' @author Guangchuang Yu
setGeneric(
    "heatplot",
    function(
        x,
        showCategory = 30,
        showTop = NULL,
        symbol = "rect",
        foldChange = NULL,
        pvalue = NULL,
        label_format = 30,
        pathway_id = NULL,
        layer = NULL,
        value = c("score", "abs_score", "share", "contribution"),
        ...
    ) {
        standardGeneric("heatplot")
    }
)

#' Volcano plot for enrichment result
#'
#'
#' @title volplot
#' @rdname volplot
#' @param x enrichment result.
#' @param color selected variable to color the dots
#' @param xintercept value to set x-intercept
#' @param yintercept value to set y-intercept
#' @param showCategory number of most significant enriched terms or selected terms to
#'     display determined by the variable selected to color the dots
#' @param label_format a numeric value setting the wrap length, alternatively a
#'     custom function to format axis labels.
#' @param ... Additional parameters
#' @export
#' @return ggplot object
#' @examples
#' \dontrun{
#' library(DOSE)
#' data(geneList)
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' volplot(x)
#' }
#' @author Guangchuang Yu
setGeneric(
    "volplot",
    function(
        x,
        color = "zScore",
        xintercept = 1,
        yintercept = 2,
        showCategory = 5,
        label_format = 30,
        ...
    ) {
        standardGeneric("volplot")
    }
)

#' Ridgeline plot for GSEA result
#'
#'
#' @title ridgeplot
#' @rdname ridgeplot
#' @param x gseaResult object
#' @param showCategory number of categories to display or a vector of terms.
#' @param fill one of "pvalue", "p.adjust", "qvalue"
#' @param core_enrichment whether to use only core_enriched genes
#' @param label_format a numeric value setting the wrap length, alternatively a
#' custom function to format axis labels.
#' @param ... additional parameters.
#' @return ggplot object
#' @export
#' @examples
#' \donttest{
#' library(DOSE)
#' data(geneList)
#' x <- gseDO(geneList)
#' ridgeplot(x)
#' }
#' @author Guangchuang Yu
setGeneric(
    "ridgeplot",
    function(
        x,
        showCategory = 30,
        fill = "p.adjust",
        core_enrichment = TRUE,
        label_format = 30,
        ...
    ) {
        standardGeneric("ridgeplot")
    }
)


#' upsetplot method generics
#'
#'
#' @docType methods
#' @name upsetplot
#' @rdname upsetplot-methods
#' @title upsetplot method
#' @param x object
#' @param n number of categories to be plotted
#' @param type one of 'boxplot' or 'violin' for `gseaResult` / `mnseaResult`
#' @param layer Optional `mnsea` layer. When `NULL`, use collapsed scores.
#' @param value score summary to display for overlapping features.
#' @param core_enrichment logical. Should only core mnsea features be used?
#' @param ... additional parameters
#' @return plot
#' @export
#' @author Guangchuang Yu
setGeneric(
    "upsetplot",
    function(
        x,
        n = 10,
        type = "boxplot",
        layer = NULL,
        value = c("score", "abs_score"),
        core_enrichment = FALSE,
        ...
    ) {
        standardGeneric("upsetplot")
    }
)


#' Functional grouping tree diagram for enrichment result of
#' over-representation test or gene set enrichment analysis.
#'
#'
#' This function visualizes gene sets as a tree.
#' Gene sets with high similarity tend to cluster together, making it easier
#' for interpretation.
#' @title treeplot
#' @rdname treeplot
#' @param x enrichment result.
#' @param showCategory number of enriched terms to display
#' @param color variable used to color enriched terms, e.g. pvalue,
#' p.adjust or qvalue
#' @param label_format a numeric value setting the wrap length, alternatively a
#' custom function to format axis labels.
#' @param ... additional parameters
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'     library(clusterProfiler)
#'     library(org.Hs.eg.db)
#'     library(enrichplot)
#'     library(GOSemSim)
#'     library(ggplot2)
#'     library(DOSE)
#'     data(geneList)
#'     gene <- names(geneList)[abs(geneList) > 2]
#'     ego <- enrichGO(gene  = gene,
#'         universe      = names(geneList),
#'         OrgDb         = org.Hs.eg.db,
#'         ont           = "BP",
#'         pAdjustMethod = "BH",
#'         pvalueCutoff  = 0.01,
#'         qvalueCutoff  = 0.05,
#'         readable      = TRUE)
#'     d <- godata('org.Hs.eg.db', ont="BP")
#'     ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
#'     treeplot(ego2, showCategory = 30)
#'     # use `hilight = FALSE` to remove ggtree::geom_hilight() layer.
#'     treeplot(ego2, showCategory = 30, hilight = FALSE)
#'     # use `offset` parameter to adjust the distance of bar and tree.
#'     treeplot(ego2, showCategory = 30, hilight = FALSE, offset = rel(1.5))
#'     # use `offset_tiplab` parameter to adjust the distance of nodes and branches.
#'     treeplot(ego2, showCategory = 30, hilight = FALSE, offset_tiplab = rel(1.5))
#'     keep <- rownames(ego2@termsim)[c(1:10, 16:20)]
#'     keep
#'     treeplot(ego2, showCategory = keep)
#'     treeplot(ego2, showCategory = 20,
#'         group_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"))
#'     # It can also graph compareClusterResult
#'     data(gcSample)
#'     xx <- compareCluster(gcSample, fun="enrichKEGG",
#'                          organism="hsa", pvalueCutoff=0.05)
#'     xx <- pairwise_termsim(xx)
#'     treeplot(xx)
#'
#'     # use `geneClusterPanel` to change the gene cluster panel.
#'     treeplot(xx, geneClusterPanel = "dotplot")
#'
#'     treeplot(xx, geneClusterPanel = "pie")
#'    }
setGeneric("treeplot", function(x, ...) {
    standardGeneric("treeplot")
})

#' Similarity Space Plot for enrichment analysis
#'
#' Creates 2D visualization of enrichment results using dimension reduction
#' techniques to show relationships between terms based on similarity.
#'
#' @title ssplot
#' @rdname ssplot
#' @inheritParams emapplot
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#'     library(clusterProfiler)
#'     library(org.Hs.eg.db)
#'     library(enrichplot)
#'     library(GOSemSim)
#'     library(DOSE)
#'     data(geneList)
#'     gene <- names(geneList)[abs(geneList) > 2]
#'     ego <- enrichGO(gene  = gene,
#'         universe      = names(geneList),
#'         OrgDb         = org.Hs.eg.db,
#'         ont           = "BP",
#'         pAdjustMethod = "BH",
#'         pvalueCutoff  = 0.01,
#'         qvalueCutoff  = 0.05,
#'         readable      = TRUE)
#'     d <- godata('org.Hs.eg.db', ont="BP")
#'     ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
#'     ssplot(ego2)
#' }
#' @author Guangchuang Yu
setGeneric("ssplot", function(x, ...) {
    standardGeneric("ssplot")
})

#' Manhattan plot for enrichment result
#'
#' @title manhattanplot
#' @rdname manhattanplot
#' @param x enrichment result.
#' @inheritParams enrichplot-common-params
#' @param ... additional parameters.
#' @return ggplot object
#' @export
setGeneric("manhattanplot", function(x, ...) {
    standardGeneric("manhattanplot")
})

#' Phase plot for enrichment-shift versus rewiring
#'
#' `phaseplot()` places each pathway on a two-dimensional plane:
#' enrichment strength (or shift) on the x-axis and rewiring score on
#' the y-axis. Point size and color encode leading-edge information and
#' mechanism class.
#'
#' @title phaseplot
#' @rdname phaseplot
#' @param x A `nseaResult` or `mnseaResult` object.
#' @param reference An optional reference result (another `nseaResult`,
#'   `mnseaResult`, or a results `data.frame`). When supplied,
#'   `phaseplot()` computes `delta_NES = NES - reference_NES` and uses it
#'   as the default x-axis.
#' @param selected_layer Optional layer name for `mnseaResult`.
#' @param reference_layer Optional reference layer name for `mnseaResult`.
#' @param x_axis One of `"delta_NES"` or `"NES"`. By default, `delta_NES`
#'   is used when `reference` is supplied, otherwise `NES`.
#' @param size_var One of `"leading_edge_size"` or `"leading_edge_overlap"`.
#' @param showCategory Number (or vector) of pathways to display.
#' @param ... Additional parameters passed to plot methods.
#' @return A ggplot object.
#' @export
setGeneric(
    "phaseplot",
    function(
        x,
        reference = NULL,
        selected_layer = NULL,
        reference_layer = NULL,
        x_axis = NULL,
        size_var = c("leading_edge_size", "leading_edge_overlap"),
        showCategory = 30,
        ...
    ) {
        standardGeneric("phaseplot")
    }
)

#' Pathway-specific rewiring plot
#'
#' @title rewireplot
#' @rdname rewireplot
#' @param x A `mnseaResult` object.
#' @param ... Additional parameters passed to plot methods.
#' @return A ggplot object.
#' @export
setGeneric("rewireplot", function(x, ...) {
    standardGeneric("rewireplot")
})

#' Multi-context mechanism consensus map
#'
#' `consensusmap()` creates a pathway by context heatmap. Tile fill
#' represents enrichment strength (`NES` or `delta_NES`) and point size
#' represents topology consistency (`rewiring_score` or
#' `leading_edge_overlap`). Mechanism class is shown as tile text and point
#' color.
#'
#' @title consensusmap
#' @rdname consensusmap
#' @param x A `mnseaResult`, `nseaResult`, or a named list of results.
#' @param fill_var One of `"NES"` or `"delta_NES"`; controls tile fill.
#' @param size_var One of `"rewiring_score"` or `"leading_edge_overlap"`;
#'   controls point size.
#' @param include_rewiring Logical; whether to draw topology-consistency points.
#' @param label One of `"mechanism_class"`, `"rewiring_score"`, or `"none"`.
#' @param reference An optional reference result used to compute `delta_NES`.
#'   For list input, the first element is used as the default reference.
#' @param ... Additional parameters passed to plot methods.
#' @return A ggplot object.
#' @export
setGeneric(
    "consensusmap",
    function(
        x,
        fill_var = c("NES", "delta_NES"),
        size_var = c("rewiring_score", "leading_edge_overlap"),
        include_rewiring = TRUE,
        label = c("mechanism_class", "rewiring_score", "none"),
        reference = NULL,
        ...
    ) {
        standardGeneric("consensusmap")
    }
)

#' Pathway state transition flow
#'
#' `mechanismflow()` draws pathway state transitions across contexts
#' (conditions or layers). Line width and point size encode flow
#' magnitude, while point shape encodes the mechanism state.
#'
#' @title mechanismflow
#' @rdname mechanismflow
#' @param x A `mnseaResult`, `nseaResult`, or a named list of results.
#' @param reference An optional reference result used to compute
#'   `delta_NES`. For list input, the first element is used as the default
#'   reference.
#' @param flow_var One of `"NES"`, `"delta_NES"`, or
#'   `"leading_edge_size"`; controls flow magnitude (line width and point
#'   size).
#' @param ... Additional parameters passed to plot methods.
#' @return A ggplot object.
#' @export
setGeneric(
    "mechanismflow",
    function(
        x,
        reference = NULL,
        flow_var = c("NES", "delta_NES", "leading_edge_size"),
        ...
    ) {
        standardGeneric("mechanismflow")
    }
)

