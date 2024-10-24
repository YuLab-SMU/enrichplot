##' @rdname goplot
##' @exportMethod goplot
setMethod("goplot", signature(x = "enrichResult"),
          function(x, showCategory = 10, color = "p.adjust",
                   layout = "sugiyama", geom="text", ...) {
              goplot.enrichResult(x, showCategory = showCategory,
                  color = color, layout = layout, geom = geom, ...)
          })

##' @rdname goplot
##' @exportMethod goplot
setMethod("goplot", signature(x = "gseaResult"),
          function(x, showCategory = 10, color = "p.adjust",
                   layout = "sugiyama", geom="text", ...) {
              goplot.enrichResult(x, showCategory = showCategory,
                  color = color, layout = layout, geom = geom, ...)
          })


##' @rdname goplot
##' @importFrom utils data
##' @import GOSemSim
##' @importFrom ggplot2 scale_fill_gradientn
##' @importFrom grid arrow
##' @importFrom grid unit
##' @importFrom rlang check_installed
##' @author Guangchuang Yu
goplot.enrichResult <- function(x, showCategory = 10, color = "p.adjust",
                                layout = "sugiyama", geom = "text", ...) {
    segment.size <- get_ggrepel_segsize()
    # has_package("AnnotationDbi")
    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    y <- y[1:n,]

    id <- y$ID[1:n]

    if (!exists(".GOSemSimEnv")) GOSemSim_initial()
    .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
    gotbl <- get("gotbl", envir=.GOSemSimEnv)
    if (inherits(x, "gseaResult")) {
        GOANCESTOR <- getAncestors(x@setType)
    } else {
        GOANCESTOR <- getAncestors(x@ontology)
    }
    
    anc <- GOANCESTOR[id] 
    ca <- anc[[1]]
    for (i in 2:length(anc)) {
        ca <- intersect(ca, anc[[i]])
    }

    uanc <- unique(unlist(anc))
    uanc <- uanc[!uanc %in% ca]
    dag <- gotbl[gotbl$go_id %in% unique(c(id, uanc)),]


    edge <- dag[, c(5, 1, 4)]
    node <- unique(gotbl[gotbl$go_id %in% unique(c(edge[,1], edge[,2])), 1:3])
    node$color <- x[node$go_id, color]
    node$size <- sapply(geneSets[node$go_id], length)

    g <- graph_from_data_frame(edge, directed=TRUE, vertices=node)
    E(g)$relationship <- edge[,3]

    p <- ggplot(g, layout = layout) +
        geom_edge(aes(linetype = .data$relationship),
            arrow = arrow(length = unit(2, 'mm')),
            colour="darkgrey") 

    if (geom == "label") {
        p <- p + geom_label_repel(aes(label=.data$label), segment.size = segment.size) +
            set_enrichplot_color(type = "fill", name = color, na.value="white")
    } else {
        p <- p + geom_text_repel(aes(label=.data$label), segment.size = segment.size, bg.color="white", bg.r=.1) 
    }        
        
    return(p)
}

##' @importFrom utils getFromNamespace
GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
getAncestors <- getFromNamespace("getAncestors", "GOSemSim")
