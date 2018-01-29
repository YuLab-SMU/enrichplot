## ##' show method for \code{enrichResult} instance
## ##'
## ##' @name show
## ##' @docType methods
## ##' @rdname show-methods
## ##'
## ##' @title show method
## ##' @param object A \code{enrichResult} instance.
## ##' @return message
## ##' @importFrom utils str
## ##' @importFrom methods show
## ##' @exportMethod show
## ##' @usage show(object)
## ##' @examples
## ##' library(DOSE)
## ##' data(geneList)
## ##' de <- names(geneList)[1:100]
## ##' x <- enrichDO(de)
## ##' print(x)
## ##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
## setMethod("show", signature(object="enrichResult"),
##           function (object){
##               cat("#\n# over-representation test\n#\n")
##               cat("#...@organism", "\t", object@organism, "\n")
##               cat("#...@ontology", "\t", object@ontology, "\n")
##               kt <- object@keytype
##               if (kt != "UNKNOWN") {
##                   cat("#...@keytype", "\t", kt, "\n")
##               }
##               cat("#...@gene", "\t")
##               str(object@gene)
##               cat("#...pvalues adjusted by", paste0("'", object@pAdjustMethod, "'"),
##                   paste0("with cutoff <", object@pvalueCutoff), "\n")
##               cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
##               str(object@result)
##               cat("#...Citation\n")
##               if (object@ontology == "DO" || object@ontology == "DOLite" || object@ontology == "NCG") {
##                   citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an",
##                                     "  R/Bioconductor package for Disease Ontology Semantic and Enrichment",
##                                     "  analysis. Bioinformatics 2015, 31(4):608-609", sep="\n", collapse="\n")
##               } else if (object@ontology == "Reactome") {
##                   citation_msg <- paste("  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for",
##                                         "  reactome pathway analysis and visualization. Molecular BioSystems",
##                                         "  2016, 12(2):477-479", sep="\n", collapse="\n")
##               } else {
##                   citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.",
##                                         "  clusterProfiler: an R package for comparing biological themes among",
##                                         "  gene clusters. OMICS: A Journal of Integrative Biology",
##                                         "  2012, 16(5):284-287", sep="\n", collapse="\n")
##               }
##               cat(citation_msg, "\n\n")
##           })

## ##' show method for \code{gseaResult} instance
## ##'
## ##' @name show
## ##' @docType methods
## ##' @rdname show-methods
## ##'
## ##' @title show method
## ##' @return message
## ##' @importFrom methods show
## ##' @exportMethod show
## ##' @usage show(object)
## ##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
## setMethod("show", signature(object="gseaResult"),
##           function (object){
##               params <- object@params
##               cat("#\n# Gene Set Enrichment Analysis\n#\n")
##               cat("#...@organism", "\t", object@organism, "\n")
##               cat("#...@setType", "\t", object@setType, "\n")
##               kt <- object@keytype
##               if (kt != "UNKNOWN") {
##                   cat("#...@keytype", "\t", kt, "\n")
##               }

##               cat("#...@geneList", "\t")
##               str(object@geneList)
##               cat("#...nPerm", "\t", params$nPerm, "\n")
##               cat("#...pvalues adjusted by", paste0("'", params$pAdjustMethod, "'"),
##                   paste0("with cutoff <", params$pvalueCutoff), "\n")
##               cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
##               str(object@result)
##               cat("#...Citation\n")
##               if (object@setType == "DO" || object@setType == "DOLite" || object@setType == "NCG") {
##                   citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an",
##                                     "  R/Bioconductor package for Disease Ontology Semantic and Enrichment",
##                                     "  analysis. Bioinformatics 2015, 31(4):608-609", sep="\n", collapse="\n")
##               } else if (object@setType == "Reactome") {
##                   citation_msg <- paste("  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for",
##                                         "  reactome pathway analysis and visualization. Molecular BioSystems",
##                                         "  2016, 12(2):477-479", sep="\n", collapse="\n")
##               } else {
##                   citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.",
##                                         "  clusterProfiler: an R package for comparing biological themes among",
##                                         "  gene clusters. OMICS: A Journal of Integrative Biology",
##                                         "  2012, 16(5):284-287", sep="\n", collapse="\n")
##               }
##               cat(citation_msg, "\n\n")
##           }
##           )

