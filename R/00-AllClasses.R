##' Class "enrichResult"
##' This class represents the result of enrichment analysis.
##'
##'
##' @name enrichResult-class
##' @aliases enrichResult-class
##'   show,enrichResult-method
##'
##' @docType class
##' @slot result enrichment analysis
##' @slot pvalueCutoff pvalueCutoff
##' @slot pAdjustMethod pvalue adjust method
##' @slot qvalueCutoff qvalueCutoff
##' @slot organism only "human" supported
##' @slot ontology biological ontology
##' @slot gene Gene IDs
##' @slot keytype Gene ID type
##' @slot universe background gene
##' @slot gene2Symbol mapping gene to Symbol
##' @slot geneSets gene sets
##' @slot readable logical flag of gene ID in symbol or not.
##' @exportClass enrichResult
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @keywords classes
setClass("enrichResult",
         representation=representation(
             result         = "data.frame",
             pvalueCutoff   = "numeric",
             pAdjustMethod  = "character",
             qvalueCutoff   = "numeric",
             organism       = "character",
             ontology       = "character",
             gene           = "character",
             keytype        = "character",
             universe       = "character",
             gene2Symbol    = "character",
             geneSets       = "list",
             readable       = "logical"
             ),
         prototype=prototype(readable = FALSE)
         )



##' Class "gseaResult"
##' This class represents the result of GSEA analysis
##'
##'
##' @name gseaResult-class
##' @aliases gseahResult-class
##'   show,gseaResult-method
##'
##' @docType class
##' @slot result GSEA anaysis
##' @slot organism organism
##' @slot setType setType
##' @slot geneSets geneSets
##' @slot geneList order rank geneList
##' @slot keytype ID type of gene
##' @slot permScores permutation scores
##' @slot params parameters
##' @slot gene2Symbol gene ID to Symbol
##' @slot readable whether convert gene ID to symbol
##' @exportClass gseaResult
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{gseaplot}}
##' @keywords classes
setClass("gseaResult",
         representation   = representation(
             result          = "data.frame",
             organism        = "character",
             setType         = "character",
             geneSets        = "list",
             geneList        = "numeric",
             keytype         = "character",
             permScores      = "matrix",
             params          = "list",
             gene2Symbol     = "character",
             readable        = "logical"
         )
         )

