## Importers for enrichment-analysis results produced by external tools
## (enrichr, g:Profiler, WebGestalt, fgsea).  Each importer maps the
## tool-specific columns to the canonical schema and delegates object
## construction to the 'enrichit' constructors, which are re-exported
## in R/reexport.R.

##' Import enrichr results
##'
##' Converts the output of \code{enrichR::enrichr()} (a list of result
##' tables, one per database, or a single result table) into an
##' \code{enrichResult} object compatible with 'enrichplot' visualizations.
##' The \code{gene} vector used in the analysis should be supplied for
##' exact \code{GeneRatio} values; \code{universe} enables
##' \code{BgRatio}/\code{RichFactor} calculation.
##' @title import_enrichr
##' @param x result table (data.frame) or the list of tables returned by
##' \code{enrichR::enrichr()}
##' @param db if \code{x} is a list, the name or index of the database
##' table to import
##' @param gene query gene vector used in the analysis
##' @param universe background gene vector used in the analysis
##' @param ontology,organism,keytype metadata stored in the object slots
##' @param pAdjustMethod method used to compute \code{p.adjust} when the
##' table has no adjusted p-value column
##' @param ... additional parameters passed to
##' \code{\link[enrichit]{as_enrichResult}}
##' @return An \code{enrichResult} object
##' @export
import_enrichr <- function(
    x,
    db = 1,
    gene = NULL,
    universe = NULL,
    ontology = "UNKNOWN",
    organism = "UNKNOWN",
    keytype = "UNKNOWN",
    pAdjustMethod = "BH",
    ...
) {
    df <- .pick_enrichr_df(x, db)
    df <- as.data.frame(df)
    .check_columns(df, c("Term", "P.value", "Adjusted.P.value", "Genes"), "enrichr")

    term <- as.character(df$Term)
    m_go <- regexpr("GO:[0-9]{7}", term)
    id <- term
    id[m_go > 0] <- regmatches(term, m_go)
    description <- trimws(gsub("GO:[0-9]{7}", "", term))
    description <- gsub("^[;:()\\s-]+|[;:()\\s-]+$", "", description, perl = TRUE)
    description[!nzchar(description)] <- term[!nzchar(description)]

    geneID <- gsub(";\\s*", "/", as.character(df$Genes))
    out <- data.frame(
        ID = id,
        Description = description,
        pvalue = as.numeric(df$P.value),
        p.adjust = as.numeric(df$Adjusted.P.value),
        geneID = geneID,
        stringsAsFactors = FALSE
    )

    if (!is.null(df$Overlap)) {
        k <- suppressWarnings(as.integer(sub("/.*$", "", as.character(df$Overlap))))
        m <- suppressWarnings(as.integer(sub("^.*/", "", as.character(df$Overlap))))
        if (!all(is.na(k))) {
            out$Count <- k
            if (!is.null(universe) && all(!is.na(m))) {
                out$BgRatio <- paste0(m, "/", length(universe))
            }
        }
    }

    enrichit::as_enrichResult(
        out,
        gene = gene,
        universe = universe,
        ontology = ontology,
        organism = organism,
        keytype = keytype,
        pAdjustMethod = pAdjustMethod,
        ...
    )
}

##' Import g:Profiler (gprofiler2) results
##'
##' Converts the result of \code{gprofiler2::gost()} (the list it returns,
##' or its \code{result} data.frame) into an \code{enrichResult} object.
##' Only the over-representation mode is mapped in this version.
##' @title import_gprofiler2
##' @param x the object returned by \code{gprofiler2::gost()}, or its
##' \code{result} component
##' @param gene query gene vector; used only when the table has no
##' \code{query_size} column
##' @param ontology,organism,keytype metadata stored in the object slots
##' @param pAdjustMethod method used to compute \code{p.adjust} when the
##' table has no adjusted p-value column
##' @param ... additional parameters passed to
##' \code{\link[enrichit]{as_enrichResult}}
##' @return An \code{enrichResult} object
##' @export
import_gprofiler2 <- function(
    x,
    gene = NULL,
    ontology = "UNKNOWN",
    organism = "UNKNOWN",
    keytype = "UNKNOWN",
    pAdjustMethod = "BH",
    ...
) {
    if (is.list(x) && !is.data.frame(x)) {
        if (is.null(x$result)) {
            stop("x does not appear to be a gost() result (no 'result' element)")
        }
        x <- x$result
    }
    df <- as.data.frame(x)
    .check_columns(
        df,
        c("term_id", "term_name", "p_value", "intersection_size"),
        "g:Profiler"
    )

    out <- data.frame(
        ID = as.character(df$term_id),
        Description = as.character(df$term_name),
        pvalue = as.numeric(df$p_value),
        geneID = .unlist_column(df$intersections),
        Count = as.numeric(df$intersection_size),
        stringsAsFactors = FALSE
    )
    if (!is.null(df$adjusted_p_value)) {
        out$p.adjust <- as.numeric(df$adjusted_p_value)
    }
    if (!is.null(df$term_size) && !is.null(df$effective_domain_size)) {
        out$BgRatio <- paste0(
            as.integer(df$term_size), "/", as.integer(df$effective_domain_size)
        )
    }
    if (!is.null(df$query_size)) {
        out$GeneRatio <- paste0(
            as.integer(df$intersection_size), "/", as.integer(df$query_size)
        )
    }

    enrichit::as_enrichResult(
        out,
        gene = gene,
        ontology = ontology,
        organism = organism,
        keytype = keytype,
        pAdjustMethod = pAdjustMethod,
        ...
    )
}

##' Import WebGestaltR results
##'
##' Converts the result table returned by \code{WebGestaltR()} into an
##' \code{enrichResult} object.  The \code{gene} vector used in the
##' analysis should be supplied for exact \code{GeneRatio} values.
##' @title import_webgestalt
##' @param x result table (data.frame) returned by \code{WebGestaltR()}
##' @param gene query gene vector used in the analysis
##' @param universe background gene vector used in the analysis
##' @param ontology,organism,keytype metadata stored in the object slots
##' @param pAdjustMethod method used to compute \code{p.adjust} when the
##' table has no adjusted p-value column
##' @param ... additional parameters passed to
##' \code{\link[enrichit]{as_enrichResult}}
##' @return An \code{enrichResult} object
##' @export
import_webgestalt <- function(
    x,
    gene = NULL,
    universe = NULL,
    ontology = "UNKNOWN",
    organism = "UNKNOWN",
    keytype = "UNKNOWN",
    pAdjustMethod = "BH",
    ...
) {
    df <- as.data.frame(x)
    pcol <- intersect(c("rawPValue", "pValue", "PValue", "p_value"), names(df))[1]
    qcol <- intersect(c("adjPValue", "FDR", "p_adj"), names(df))[1]
    required <- c("geneSet", "description", pcol, "overlap", "userIds")
    .check_columns(df, required[!is.na(required)], "WebGestalt")

    out <- data.frame(
        ID = as.character(df$geneSet),
        Description = as.character(df$description),
        pvalue = as.numeric(df[[pcol]]),
        geneID = gsub(";\\s*", "/", as.character(df$userIds)),
        Count = as.numeric(df$overlap),
        stringsAsFactors = FALSE
    )
    if (!is.na(qcol)) {
        out$p.adjust <- as.numeric(df[[qcol]])
    }
    if (!is.null(df$size) && !is.null(universe)) {
        out$BgRatio <- paste0(as.integer(df$size), "/", length(universe))
    }

    enrichit::as_enrichResult(
        out,
        gene = gene,
        universe = universe,
        ontology = ontology,
        organism = organism,
        keytype = keytype,
        pAdjustMethod = pAdjustMethod,
        ...
    )
}

##' Import fgsea results
##'
##' Converts the result table returned by \code{fgsea::fgsea()} into a
##' \code{gseaResult} object.  The ranked statistics vector (\code{stats})
##' used as \code{fgsea()} input is required, as most GSEA visualizations
##' consume the ranked list.  Passing the \code{pathways} list used in the
##' analysis enables exact running-score plots; without it, gene sets are
##' rebuilt from the leading edge and plots are approximate.
##' @title import_fgsea
##' @param x result table returned by \code{fgsea::fgsea()}
##' @param stats named numeric vector of ranked statistics, sorted in
##' descending order (sorted automatically with a warning if not)
##' @param geneSets the \code{pathways} list used as \code{fgsea()} input
##' @param setType,organism,keytype metadata stored in the object slots
##' @param exponent,scoreType parameters used to recompute missing
##' \code{rank}/\code{leading_edge}/\code{core_enrichment} columns
##' @param pAdjustMethod method used to compute \code{p.adjust} when the
##' table has no adjusted p-value column
##' @param ... additional parameters passed to
##' \code{\link[enrichit]{as_gseaResult}}
##' @return A \code{gseaResult} object
##' @export
import_fgsea <- function(
    x,
    stats,
    geneSets = NULL,
    setType = "UNKNOWN",
    organism = "UNKNOWN",
    keytype = "UNKNOWN",
    exponent = 1,
    scoreType = "std",
    pAdjustMethod = "BH",
    ...
) {
    df <- as.data.frame(x)
    .check_columns(
        df,
        c("pathway", "ES", "pval"),
        "fgsea"
    )

    enrichit::as_gseaResult(
        df,
        geneList = stats,
        geneSets = geneSets,
        setType = setType,
        organism = organism,
        keytype = keytype,
        exponent = exponent,
        scoreType = scoreType,
        pAdjustMethod = pAdjustMethod,
        ...
    )
}

## ---- internal helpers -------------------------------------------------

.pick_enrichr_df <- function(x, db) {
    if (is.data.frame(x)) {
        return(x)
    }
    if (is.list(x)) {
        if (is.character(db)) {
            if (!db %in% names(x)) {
                stop(
                    "db '", db, "' not found; available: ",
                    paste(utils::head(names(x), 10), collapse = ", ")
                )
            }
            return(x[[db]])
        }
        if (length(db) == 1 && is.numeric(db)) {
            if (db > length(x)) {
                stop("db index out of range; x contains ", length(x), " table(s)")
            }
            return(x[[db]])
        }
    }
    stop("x must be a data.frame or the list returned by enrichR::enrichr()")
}

.check_columns <- function(df, required, what) {
    miss <- setdiff(required, names(df))
    if (length(miss) > 0) {
        stop(
            what, " output is missing required column(s): ",
            paste(miss, collapse = ", ")
        )
    }
}

.unlist_column <- function(col) {
    if (is.null(col)) {
        return(NULL)
    }
    if (is.list(col)) {
        vapply(col, function(g) paste0(unique(as.character(g)), collapse = "/"), character(1))
    } else {
        gsub("[;,]\\s*", "/", as.character(col))
    }
}
