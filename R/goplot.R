#' @rdname goplot
#' @exportMethod goplot
setMethod("goplot", signature(x = "enrichResult"),
          function(x, showCategory = 10, color = "p.adjust",
                   layout = igraph::layout_with_sugiyama, geom="text", ...) {
              goplot.enrichResult(x, showCategory = showCategory,
                  color = color, layout = layout, geom = geom, ...)
          })

#' @rdname goplot
#' @exportMethod goplot
setMethod("goplot", signature(x = "gseaResult"),
          function(x, showCategory = 10, color = "p.adjust",
                   layout = igraph::layout_with_sugiyama, geom="text", ...) {
              goplot.enrichResult(x, showCategory = showCategory,
                  color = color, layout = layout, geom = geom, ...)
          })



#' @importFrom utils data
#' @import GOSemSim
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom grid arrow
#' @importFrom grid unit
#' @importFrom rlang check_installed
#' @importFrom yulab.utils get_cache_item
goplot.enrichResult <- function(x, showCategory = 10, color = "p.adjust",
                                layout = igraph::layout_with_sugiyama, geom = "text", 
                                ID = "Description", ...) {
    segment.size <- get_ggrepel_segsize()
    require_suggested('glue', 'for `goplot()`.')
    # has_package("AnnotationDbi")
    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    y <- y[1:n,]

    id <- y$ID[1:n]

    GOSemSim_initial()
    .GOSemSimEnv <- get_cache_item(".GOSemSimEnv")
    gotbl <- .GOSemSimEnv$gotbl

    if (inherits(gotbl, "character")) {
        utils::data("gotbl", package = "GOSemSim", envir = environment())
        gotbl <- get("gotbl")
    }

    if (inherits(x, "gseaResult")) {
        onto <- x@setType
    } else {
        onto <- x@ontology
    }

    if (!toupper(onto) %in% c("MF", "CC", "BP")) {
        stop("Ontology should be one of 'MF', 'CC' or 'BP'")
    }

    GOANCESTOR <- getAncestors(onto)

    anc <- GOANCESTOR[id] 
    ca <- anc[[1]]
    for (i in 2:length(anc)) {
        ca <- intersect(ca, anc[[i]])
    }

    uanc <- unique(unlist(anc))
    uanc <- uanc[!uanc %in% ca]
    dag_ids <- unique(c(id, uanc))
    dag <- gotbl[gotbl$go_id %in% dag_ids,]

    edge <- dag[, c(5, 1, 4)]
    node <- unique(gotbl[gotbl$go_id %in% dag_ids, 1:3])
    valid_vertices <- node$go_id
    edge <- edge[
        edge[, 1] %in% valid_vertices &
        edge[, 2] %in% valid_vertices,
        ,
        drop = FALSE
    ]
    node$color <- x[node$go_id, color]
    node$size <- sapply(geneSets[node$go_id], length)

    g <- graph_from_data_frame(edge, directed=TRUE, vertices=node)
    E(g)$relationship <- edge[,3]

    require_suggested('ggarchery', 'for `goplot()`.')

    position = ggarchery::position_attractsegment(
            start_shave=.03, 
            end_shave=.03,
            type_shave="proportion"
        )
    p <- ggplot(g, layout = layout) +
        geom_edge(aes(linetype = .data$relationship),
            arrow = arrow(length = unit(2, 'mm')),
            colour="darkgrey", position=position, linewidth = 1) 

    if (ID == "Description" || ID == "ID") {
        ID <- sprintf("{%s}", ID)
    } 

    if (geom == "label") {
        p <- p + geom_label_repel(aes(label= glue::glue(ID, ID=.data[['name']], Description=.data[['Term']]), 
                        fill=.data$color, segment.size = segment.size)) +
            set_enrichplot_color(type = "fill", name = color, na.value="white")
    } else {
        p <- p + geom_point(aes(color=.data$color), size=5) +
            geom_text_repel(aes(label=glue::glue(ID, ID=.data[['name']], Description=.data[['Term']])), 
                    segment.size = segment.size, bg.color="white", bg.r=.1) +
            set_enrichplot_color(type = "color", name = color, na.value="grey")
    }        
        
    return(p)
}

#' @importFrom utils getFromNamespace
GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
getAncestors <- getFromNamespace("getAncestors", "GOSemSim")
