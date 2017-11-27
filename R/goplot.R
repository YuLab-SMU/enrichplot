##' plot induced GO DAG of significant terms
##'
##'
##' @title goplot
##' @param x enrichment result. e.g. instance of gseaResult or enrichResult
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue, p.adjust or qvalue
##' @param layout layout of the map
##' @param ... additional parameter
##' @return ggplot object
##' @importFrom utils data
##' @import GOSemSim
##' @importFrom ggplot2 scale_fill_gradientn
##' @importFrom grid arrow
##' @importFrom grid unit
##' @importFrom ggraph circle
##' @importFrom ggraph geom_node_label
##' @export
##' @author guangchuang yu
goplot <- function(x, showCategory = 10, color = "p.adjust", layout = "sugiyama", ...) {
    id <- x$ID[1:showCategory]

    if (!exists(".GOSemSimEnv")) GOSemSim_initial()
    .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
    gotbl <- get("gotbl", envir=.GOSemSimEnv)
    ## data("gotbl", package="GOSemSim")
    ## gotbl <- get("gotbl", envir=.GlobalEnv)
    dag = gotbl[gotbl$go_id %in% id, ]

    cid <- dag$go_id
    cid = unique(dag$parent[!dag$parent %in% cid])
    if (length(cid) > 1 || cid != "all") {
        cdag = gotbl[gotbl$go_id %in% cid, ]
        dag <- unique(rbind(dag, cdag))
    }

    edge <- dag[, c(5, 1, 4)]
    node <- unique(gotbl[gotbl$go_id %in% unique(c(edge[,1], edge[,2])), 1:3])
    node$color <- x[node$go_id, color]

    g <- graph.data.frame(edge, directed=T, vertices=node)

    ggraph(g, layout=layout) +
        geom_edge_link(arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm')) +
        geom_node_point(size=5, aes_(color=~color)) +
        geom_node_label(aes_(label=~Term, fill=~color), repel=T) +
        theme_void() +
        scale_color_gradientn(name = color, colors=heatmap_palette, guide=guide_colorbar(reverse=TRUE)) +
        scale_fill_gradientn(name = color, colors=heatmap_palette, guide=guide_colorbar(reverse=TRUE), na.value='white')
}

GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
