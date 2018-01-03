##' @method geneInCategory enrichResult
##' @export
##' @importFrom stats setNames
geneInCategory.enrichResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))

##' @method geneInCategory gseaResult
##' @export
geneInCategory.gseaResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))

##' @method geneID enrichResult
##' @export
geneID.enrichResult <- function(x) as.character(x@result$geneID)

##' @method geneID gseaResult
##' @export
geneID.gseaResult <- function(x) as.character(x@result$core_enrichment)
