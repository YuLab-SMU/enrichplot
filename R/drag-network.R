
#' Drag the nodes of a network to update the layout of the network
#'
#' @param p the network diagram as a ggplot/gg/ggraph object. 
#' @param g an corresponding igraph object. Default is to extract from the 'ggraph' attribute.
#'
#' @importFrom igraph tk_coords
#' @importFrom igraph tkplot
#' @importFrom igraph V
#' @export
#' @return an updated ggplot/gg/ggraph object
#' @examples
#' \dontrun{
#' library(igraph)
#' library(ggraph)
#' 
#' flow_info <- data.frame(from = c(1,2,3,3,4,5,6),
#'                         to = c(5,5,5,6,7,6,7))
#' g = graph_from_data_frame(flow_info)
#' p <- ggraph(g, layout='nicely') + geom_node_point() + geom_edge_link() 
#' pp <- drag_network(p)
#' }
drag_network <- function(p, g = NULL) {
    if (is.null(g)) g = attr(p$data, "graph")
    nodeColor <- V(g)$color
    # If the color of the node is continuous, an error will be reported.
    # So assign the color to black.
    if (length(.tkplot.convert.color(nodeColor)) == 0) {
        V(g)$color <- "black"
    }

    x <- tkplot(g)
    # if use 'igraph:::.tkplot.get(x)', the plot will upside down.
    y <- tk_coords(x)
    while(TRUE) {
        y2 <- tryCatch(tk_coords(x), error=function(e) NULL)
        if (is.null(y2)) break
        if (identical(y, y2)) {
            Sys.sleep(1)
            next
        } else {
            y <- y2
            p$data$x <- y[,1]
            p$data$y <- y[,2]
            print(p)
        }
    }
    invisible(p)
}


.tkplot.convert.color <- getFromNamespace(".tkplot.convert.color", "igraph")
