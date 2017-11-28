##' visualize analyzing result of GSEA
##'
##' plotting function for gseaResult
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_linerange
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggplotGrob
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 rel
##' @importFrom cowplot plot_grid
##' @param gseaResult gseaResult object
##' @param geneSetID geneSet ID
##' @param by one of "runningScore" or "position"
##' @param title plot title
##' @param color color of line segments
##' @param color.line color of running enrichment score line
##' @param color.vline color of vertical line which indicating the maximum/minimal running enrichment score
##' @return ggplot2 object
##' @export
##' @examples
##' library(DOSE)
##' data(geneList)
##' x <- gseDO(geneList)
##' gseaplot(x, geneSetID=1)
##' @author Yu Guangchuang
gseaplot <- function (gseaResult, geneSetID, by = "all", title = "", color='black', color.line="green", color.vline="#FA5860"){
    by <- match.arg(by, c("runningScore", "preranked", "all"))
    gsdata <- gsInfo(gseaResult, geneSetID)
    p <- ggplot(gsdata, aes_(x = ~x)) +
        theme_dose() + xlab("Position in the Ranked List of Genes")
    if (by == "runningScore" || by == "all") {
        p.res <- p + geom_linerange(aes_(ymin=~ymin, ymax=~ymax), color=color)
        p.res <- p.res + geom_line(aes_(y = ~runningScore), color=color.line, size=1)
        enrichmentScore <- gseaResult@result[geneSetID, "enrichmentScore"]
        es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
        p.res <- p.res + geom_vline(data = es.df, aes_(xintercept = ~es),
                                    colour = color.vline, linetype = "dashed")
        p.res <- p.res + ylab("Running Enrichment Score")
        p.res <- p.res + geom_hline(yintercept = 0)
    }
    if (by == "preranked" || by == "all") {
        df2 <- data.frame(x = which(p$data$position == 1))
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), color=color)
        p.pos <- p.pos + ylab("Ranked list metric") + xlim(0, length(p$data$geneList))
    }
    if (by == "runningScore")
        return(p.res + ggtitle(title))
    if (by == "preranked")
        return(p.pos + ggtitle(title))

    p.pos <- p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank())
    p.pos <- p.pos + ggtitle(title) +
        theme(plot.title=element_text(hjust=0.5, size=rel(2)))
    plot_grid(p.pos, p.res, ncol=1, align="v")
}


##' extract gsea result of selected geneSet
##'
##'
##' @title gsInfo
##' @param object gseaResult object
##' @param geneSetID gene set ID
##' @return data.frame
##' @author Guangchuang Yu
## @export
gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin=0
    df$ymax=0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList

    return(df)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")
