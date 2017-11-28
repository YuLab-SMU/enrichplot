##' plot induced GO DAG of significant terms
##'
##'
##' @title goplot
##' @param x enrichment result. e.g. instance of gseaResult or enrichResult
##' @param showCategory number of enriched terms to display
##' @param color variable that used to color enriched terms, e.g. pvalue, p.adjust or qvalue
##' @param layout layout of the map
##' @param geom label geom, one of 'label' or 'text'
##' @param ... additional parameter
##' @return ggplot object
##' @importFrom utils data
##' @import GOSemSim
##' @importFrom ggplot2 scale_fill_gradientn
##' @importFrom grid arrow
##' @importFrom grid unit
##' @importFrom ggraph circle
##' @importFrom ggraph geom_node_label
##' @importFrom AnnotationDbi mget
##' @export
##' @author guangchuang yu
goplot <- function(x, showCategory = 10, color = "p.adjust", layout = "sugiyama", geom = "text", ...) {
    if (!inherits(x, "gseaResult") && !inherits(x, "enrichResult"))
        stop("object not supported...")

    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    y <- y[1:n,]

    id <- y$ID[1:n]

    if (!exists(".GOSemSimEnv")) GOSemSim_initial()
    .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
    gotbl <- get("gotbl", envir=.GOSemSimEnv)
    ## data("gotbl", package="GOSemSim")
    ## gotbl <- get("gotbl", envir=.GlobalEnv)

    ## dag = gotbl[gotbl$go_id %in% id, ]

    ## cid <- dag$go_id
    ## cid = unique(dag$parent[!dag$parent %in% cid])
    ## if (length(cid) > 1 || cid != "all") {
    ##     cdag = gotbl[gotbl$go_id %in% cid, ]
    ##     dag <- unique(rbind(dag, cdag))
    ## }

    GOANCESTOR <- getAncestors(x@ontology)
    anc <- AnnotationDbi::mget(id, GOANCESTOR)
    ca <- anc[[1]]
    for (i in 2:length(anc)) {
        ca <- intersect(ca, anc[[i]])
    }
    ## mrca <- ca[1]
    uanc <- unique(unlist(anc))
    uanc <- uanc[!uanc %in% ca]
    ## uanc <- uanc[sapply(uanc, function(g) sum(sapply(anc, function(y) g %in% y))) > 1]

    dag <- gotbl[gotbl$go_id %in% unique(c(id, uanc)),]


    edge <- dag[, c(5, 1, 4)]
    node <- unique(gotbl[gotbl$go_id %in% unique(c(edge[,1], edge[,2])), 1:3])
    node$color <- x[node$go_id, color]
    node$size <- sapply(geneSets[node$go_id], length)

    g <- graph.data.frame(edge, directed=TRUE, vertices=node)
    E(g)$relationship <- edge[,3]

    p <- ggraph(g, layout=layout) +
        geom_edge_link(aes_(linetype = ~relationship), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm'), colour="darkgrey") +
        geom_node_point(size = 5, aes_(color=~color)) +
        theme_void() +
        scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE))

    if (geom == "label") {
        p <- p + geom_node_label(aes_(label=~Term, fill=~color), repel=TRUE) +
            scale_fill_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE), na.value='white')
    } else {
        p <- p + geom_node_text(aes_(label=~Term), repel=TRUE)
    }
    return(p)
}

GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
getAncestors <- getFromNamespace("getAncestors", "GOSemSim")
