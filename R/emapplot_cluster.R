##' Functional grouping network diagram for enrichment result of
##' over-representation test or gene set enrichment analysis
##'
##' This function has been replaced by `emapplot`.
##'
##' @param x enrichment result
##' @param ... additional parameters. Please refer to: \link{emapplot}.
##'
##' @return ggplot2 object
##' @export
emapplot_cluster <- function(x, ...){

    .Deprecated(msg = "'emapplot_cluster' has been replaced by `emapplot`, 
        it will be removed in the next version")
    emapplot(x = x,
            group_category = TRUE,
            node_label  = "group",
            with_edge = TRUE, ...)
}















