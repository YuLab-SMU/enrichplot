## Classes used for results imported from external enrichment tools.
##
## enrichit's default show() method prints the clusterProfiler citation for
## every enrichResult/gseaResult. Imported results did not come from that
## workflow, so keep the normal class inheritance for plotting while using a
## citation-neutral display method.

methods::setClass("externalEnrichResult", contains = "enrichResult")
methods::setClass("externalGseaResult", contains = "gseaResult")

methods::setMethod(
    "show",
    "externalEnrichResult",
    function(object) {
        cat("#\n# imported over-representation result\n#\n")
        cat("#...@organism\t", object@organism, "\n", sep = "")
        cat("#...@ontology\t", object@ontology, "\n", sep = "")
        cat("#...", nrow(object@result), "imported terms found\n", sep = "")
        if (nrow(object@result) > 0) {
            utils::str(object@result)
        }
    }
)

methods::setMethod(
    "show",
    "externalGseaResult",
    function(object) {
        cat("#\n# imported gene-set enrichment result\n#\n")
        cat("#...@organism\t", object@organism, "\n", sep = "")
        cat("#...@setType\t", object@setType, "\n", sep = "")
        cat("#...", nrow(object@result), "imported gene sets found\n", sep = "")
        if (nrow(object@result) > 0) {
            utils::str(object@result)
        }
    }
)

.as_external_enrich_result <- function(object) {
    methods::as(object, "externalEnrichResult")
}

.as_external_gsea_result <- function(object) {
    methods::as(object, "externalGseaResult")
}
