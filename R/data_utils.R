#' Data processing utility functions for enrichplot package
#'
#' This file contains data manipulation and processing helper functions

#' Update showCategory parameter
#'
#' @param x input object
#' @param showCategory category specification
#' @return updated category specification
#' @noRd
update_n <- function(x, showCategory) {
    ## Input validation
    check_input(x, arg_name = "x")
    check_input(showCategory, arg_name = "showCategory")

    if (!is.numeric(showCategory)) {
        if (inherits(x, 'list')) {
            showCategory <- showCategory[showCategory %in% names(x)]
        } else {
            if (!all(c("ID", "Description") %in% colnames(as.data.frame(x)))) {
                yulab_abort(
                    "Input data must have 'ID' and 'Description' columns",
                    class = "missing_column_error"
                )
            }
            mapping <- get_term_mapping(x)
            valid_terms <- unique(c(
                mapping$ID,
                mapping$Description,
                mapping$label
            ))
            showCategory <- showCategory[showCategory %in% valid_terms]
        }
        return(showCategory)
    }

    n <- showCategory
    if (inherits(x, 'list')) {
        nn <- length(x)
    } else {
        nn <- nrow(x)
    }
    if (nn < n) {
        yulab_warn(
            paste0(
                "showCategory (",
                n,
                ") is larger than available items (",
                nn,
                "). Using ",
                nn
            ),
            class = "showCategory_warning"
        )
        n <- nn
    }

    return(n)
}

get_term_mapping <- function(x) {
    y <- as.data.frame(x)
    if (!all(c("ID", "Description") %in% colnames(y))) {
        yulab_abort(
            "Input data must have 'ID' and 'Description' columns",
            class = "missing_column_error"
        )
    }

    ids <- as.character(y$ID)
    descriptions <- as.character(y$Description)
    labels <- descriptions
    duplicated_desc <- duplicated(descriptions) |
        duplicated(descriptions, fromLast = TRUE)
    labels[duplicated_desc] <- paste0(
        descriptions[duplicated_desc],
        " [",
        ids[duplicated_desc],
        "]"
    )

    data.frame(
        ID = ids,
        Description = descriptions,
        label = labels,
        stringsAsFactors = FALSE
    )
}

resolve_term_rows <- function(x, terms) {
    mapping <- get_term_mapping(x)
    idx <- integer()

    for (term in as.character(terms)) {
        id_hits <- which(mapping$ID == term)
        if (length(id_hits) > 0) {
            idx <- c(idx, id_hits)
            next
        }

        label_hits <- which(mapping$label == term)
        if (length(label_hits) > 0) {
            idx <- c(idx, label_hits)
            next
        }

        description_hits <- which(mapping$Description == term)
        if (length(description_hits) > 0) {
            idx <- c(idx, description_hits)
        }
    }

    unique(idx)
}

get_term_labels <- function(x, ids) {
    mapping <- get_term_mapping(x)
    idx <- match(as.character(ids), mapping$ID)
    labels <- mapping$label[idx]
    names(labels) <- as.character(ids)
    labels
}

get_geneSet_labels <- function(geneSets) {
    labels <- attr(geneSets, "term_labels", exact = TRUE)
    if (is.null(labels)) {
        labels <- setNames(names(geneSets), names(geneSets))
    }
    labels
}

set_geneSet_labels <- function(geneSets, labels) {
    attr(geneSets, "term_labels") <- labels
    geneSets
}

attach_result_attributes <- function(geneSets, result_df) {
    if (is.null(result_df) || nrow(result_df) == 0) {
        return(geneSets)
    }
    for (col in colnames(result_df)) {
        if (is.numeric(result_df[[col]])) {
            attr(geneSets, col) <- result_df[[col]]
        }
    }
    geneSets
}

as_plot_geneSets <- function(geneSets) {
    plot_geneSets <- geneSets
    names(plot_geneSets) <- unname(get_geneSet_labels(geneSets))
    plot_geneSets
}

select_terms <- function(x, showCategory) {
    selection <- update_n(x, showCategory)

    if (inherits(x, "list")) {
        if (is.numeric(selection)) {
            geneSets <- x[seq_len(selection)]
        } else {
            geneSets <- x[selection]
        }
        labels <- get_geneSet_labels(geneSets)
        return(list(
            selection = selection,
            result = NULL,
            geneSets = geneSets,
            ids = names(geneSets),
            labels = labels
        ))
    }

    y <- as.data.frame(x)
    if (is.numeric(selection)) {
        selected <- y[seq_len(selection), , drop = FALSE]
    } else {
        idx <- resolve_term_rows(x, selection)
        selected <- y[idx, , drop = FALSE]
    }

    labels <- get_term_labels(x, selected$ID)
    geneSets <- geneInCategory(x)[as.character(selected$ID)]
    geneSets <- set_geneSet_labels(geneSets, labels)

    list(
        selection = selection,
        result = selected,
        geneSets = geneSets,
        ids = as.character(selected$ID),
        labels = labels
    )
}

#' Extract gene sets from enrichment result
#'
#' @param x enrichment result object
#' @param n number of categories or specific categories
#' @return gene sets list
#' @noRd
extract_geneSets <- function(x, n) {
    select_terms(x, n)$geneSets
}

#' Make fold change data readable
#'
#' @param x enrichment result object
#' @param foldChange fold change vector
#' @return readable fold change vector
#' @noRd
fc_readable <- function(x, foldChange = NULL) {
    if (is.null(foldChange)) {
        return(NULL)
    }

    if (x@readable && x@keytype != "SYMBOL") {
        gid <- names(foldChange)
        if (is(x, 'gseaResult')) {
            ii <- gid %in% names(x@geneList)
        } else {
            ii <- gid %in% x@gene
        }
        gid[ii] <- x@gene2Symbol[gid[ii]]
        names(foldChange) <- gid
    }
    return(foldChange)
}

#' Calculate overlap ratio between two gene sets
#'
#' @param x first gene set
#' @param y second gene set
#' @return Jaccard similarity coefficient
#' @noRd
overlap_ratio <- function(x, y) {
    x <- unique(unlist(x))
    y <- unique(unlist(y))
    length(intersect(x, y)) / length(union(x, y))
}

#' Calculate Jaccard similarity matrix
#'
#' @param gsetlist list of gene sets
#' @param id gene set IDs
#' @param name gene set names
#' @return similarity matrix
#' @noRd
.cal_jc_similarity <- function(gsetlist, id = NULL, name = NULL) {
    if (is.null(id)) {
        id <- names(gsetlist)
    }
    n <- length(id)
    w <- matrix(NA, nrow = n, ncol = n)
    if (is.null(name)) {
        name <- id
    }
    colnames(w) <- rownames(w) <- name

    # Vectorized computation: precompute all gene sets
    gsets <- lapply(gsetlist[id], unique)

    # Use outer function for vectorized computation
    jc_matrix <- outer(seq_len(n), seq_len(n), Vectorize(function(i, j) {
        if (i == j) {
            return(1)
        }
        overlap_ratio(gsets[[i]], gsets[[j]])
    }))

    # Ensure symmetry
    jc_matrix[lower.tri(jc_matrix)] <- t(jc_matrix)[lower.tri(t(jc_matrix))]

    colnames(jc_matrix) <- rownames(jc_matrix) <- name

    return(jc_matrix)
}

#' Prepare pie data for genes in cnetplot (compareClusterResult only)
#'
#' @param y data frame from compareClusterResult
#' @return pie data
#' @importFrom rlang check_installed
#' @noRd
prepare_pie_gene <- function(y) {
    ## Input validation
    check_input(y, type = "data.frame", arg_name = "y")

    require_suggested(c('tibble', 'tidyr'), 'for `prepare_pie_gene()`.')
    gene_pie <- tibble::as_tibble(y[, c("Cluster", "Description", "geneID")])
    gene_pie$geneID <- strsplit(gene_pie$geneID, '/')
    gene_pie2 <- as.data.frame(tidyr::unnest(gene_pie, cols = geneID))
    gene_pie2 <- unique(gene_pie2)
    prepare_pie_data(gene_pie2, pie = "equal", type = "gene")
}

#' Prepare pie data for categories in cnetplot/emapplot
#'
#' @param enrichDf enrichment data frame
#' @param pie proportion type (equal, count, Count)
#' @return pie data matrix
#' @noRd
prepare_pie_category <- function(enrichDf, pie = "equal") {
    pie <- match.arg(pie, c("equal", "count", "Count"))
    if (pie == "count") {
        pie <- "Count"
    }

    pie_data <- enrichDf[, c("Cluster", "Description", "Count")]
    pie_data[, "Description"] <- as.character(pie_data[, "Description"])
    prepare_pie_data(pie_data, pie = pie)
}

#' Prepare pie data matrix
#'
#' @param pie_data input data
#' @param pie proportion type
#' @param type gene or category
#' @return pie data matrix
#' @noRd
prepare_pie_data <- function(pie_data, pie = "equal", type = "category") {
    if (type == "category") {
        ID_unique <- unique(pie_data[, 2])
    } else {
        ID_unique <- unique(pie_data[, 3])
    }

    Cluster_unique <- unique(pie_data[, 1])
    ID_Cluster_mat <- matrix(
        0,
        nrow = length(ID_unique),
        ncol = length(Cluster_unique)
    )
    rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringsAsFactors = FALSE)

    if (pie == "Count") {
        # Vectorized matrix indexing
        idx <- cbind(
            match(pie_data[, 2], rownames(ID_Cluster_mat)),
            match(pie_data[, 1], colnames(ID_Cluster_mat))
        )
        ID_Cluster_mat[idx] <- pie_data[, 3]
        # Convert all columns to numeric at once
        ID_Cluster_mat[] <- lapply(ID_Cluster_mat, as.numeric)
        return(ID_Cluster_mat)
    }

    # Vectorized matrix indexing for equal pie
    if (type == "category") {
        idx <- cbind(
            match(pie_data[, 2], rownames(ID_Cluster_mat)),
            match(pie_data[, 1], colnames(ID_Cluster_mat))
        )
    } else {
        idx <- cbind(
            match(pie_data[, 3], rownames(ID_Cluster_mat)),
            match(pie_data[, 1], colnames(ID_Cluster_mat))
        )
    }
    ID_Cluster_mat[idx] <- 1
    return(ID_Cluster_mat)
}

#' Convert compareClusterResult to data frame
#'
#' @param x compareClusterResult object
#' @param ... additional parameters
#' @return data frame
#' @export
#' @method as.data.frame compareClusterResult
as.data.frame.compareClusterResult <- function(x, ...) {
    as.data.frame(x@compareClusterResult, ...)
}



