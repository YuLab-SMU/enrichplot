#' Compute a transparent rewiring score for network-enrichment results
#'
#' The first version uses `1 - leading-edge Jaccard overlap` as the unified
#' rewiring score. A single `nseaResult` without a reference cannot produce a
#' meaningful rewiring score and returns `NA_real_`.
#'
#' @param x A `nseaResult`, `mnseaResult`, or a named list of such results.
#' @param reference_layer Optional reference layer for `mnseaResult`.
#' @param selected_layer Optional layer to compare against the reference.
#' @param threshold_abs_score Optionally threshold used to call a feature as
#'   contributing to the leading edge.
#' @param ... Additional arguments reserved for future composite metrics.
#' @return A data.frame with at least `ID`, `Description`, `layer`,
#'   `rewiring_score` and `leading_edge_overlap`.
#' @noRd
compute_rewiring_score <- function(
    x,
    reference = NULL,
    reference_layer = NULL,
    selected_layer = NULL,
    threshold_abs_score = 0,
    ...
) {
    if (!is.null(reference)) {
        if (!inherits(reference, c("nseaResult", "mnseaResult"))) {
            stop("`reference` must be a nseaResult or mnseaResult when computing cross-object rewiring.")
        }
        x_sets <- .result_feature_sets(x, layer = selected_layer)
        ref_sets <- .result_feature_sets(reference, layer = reference_layer)
        ids <- intersect(names(x_sets), names(ref_sets))
        if (length(ids) == 0) {
            return(data.frame(
                ID = character(0),
                Description = character(0),
                layer = character(0),
                leading_edge_overlap = numeric(0),
                rewiring_score = numeric(0),
                centrality_shift = numeric(0)
            ))
        }
        overlap <- vapply(ids, function(id) {
            a <- x_sets[[id]]
            b <- ref_sets[[id]]
            if (length(a) == 0 || length(b) == 0) return(0)
            length(intersect(a, b)) / length(union(a, b))
        }, numeric(1))
        desc <- get_term_labels(x, ids)
        return(data.frame(
            ID = ids,
            Description = unname(desc),
            layer = if (is.null(selected_layer)) "collapsed" else as.character(selected_layer),
            leading_edge_overlap = unname(overlap),
            rewiring_score = 1 - unname(overlap),
            centrality_shift = NA_real_,
            stringsAsFactors = FALSE
        ))
    }

    if (is.list(x) && !inherits(x, c("nseaResult", "mnseaResult"))) {
        if (is.null(names(x))) {
            stop("A list input to compute_rewiring_score() must be named.")
        }
        if (length(x) < 2) {
            stop("At least two contexts are required to compute rewiring.")
        }
        reference <- x[[1]]
        output_list <- lapply(seq_along(x), function(i) {
            context <- names(x)[i]
            res <- compute_rewiring_score(
                x[[i]],
                reference = reference,
                reference_layer = NULL,
                selected_layer = NULL,
                threshold_abs_score = threshold_abs_score
            )
            res$context <- context
            res
        })
        return(do.call(rbind, output_list))
    }

    if (inherits(x, "mnseaResult")) {
        if (is.null(reference_layer) && is.null(selected_layer)) {
            ids <- as.character(x@result$ID)
            desc <- get_term_labels(x, ids)
            return(data.frame(
                ID = ids,
                Description = unname(desc),
                layer = "collapsed",
                leading_edge_overlap = NA_real_,
                rewiring_score = NA_real_,
                centrality_shift = NA_real_,
                stringsAsFactors = FALSE
            ))
        }

        selected_feature_layer <- if (is.null(selected_layer)) "collapsed" else selected_layer
        feature_sets <- if (is.null(selected_layer)) {
            .mnsea_feature_sets(x, layer = NULL)
        } else {
            .mnsea_feature_sets(x, layer = selected_layer)
        }
        ref_sets <- if (is.null(reference_layer)) {
            .mnsea_feature_sets(x, layer = NULL)
        } else {
            .mnsea_feature_sets(x, layer = reference_layer)
        }
        ids <- intersect(names(feature_sets), names(ref_sets))
        if (length(ids) == 0) {
            return(data.frame(
                ID = character(0),
                Description = character(0),
                layer = character(0),
                leading_edge_overlap = numeric(0),
                rewiring_score = numeric(0),
                centrality_shift = numeric(0)
            ))
        }
        overlap <- vapply(ids, function(id) {
            a <- feature_sets[[id]]
            b <- ref_sets[[id]]
            if (length(a) == 0 || length(b) == 0) return(0)
            length(intersect(a, b)) / length(union(a, b))
        }, numeric(1))
        desc <- get_term_labels(x, ids)
        data.frame(
            ID = ids,
            Description = unname(desc),
            layer = selected_feature_layer,
            leading_edge_overlap = unname(overlap),
            rewiring_score = 1 - unname(overlap),
            centrality_shift = NA_real_,
            stringsAsFactors = FALSE
        )
    } else if (inherits(x, "nseaResult")) {
        ids <- as.character(x@result$ID)
        desc <- get_term_labels(x, ids)
        data.frame(
            ID = ids,
            Description = unname(desc),
            layer = "network",
            leading_edge_overlap = NA_real_,
            rewiring_score = NA_real_,
            centrality_shift = NA_real_,
            stringsAsFactors = FALSE
        )
    } else {
        stop("x must be a nseaResult, mnseaResult, or a named list of them.")
    }
}

#' Extract feature-level rewiring information
#'
#' @param x A `mnseaResult`.
#' @param pathway_id Pathway ID to compare.
#' @param selected_layer Layer to compare.
#' @param reference_layer Reference layer (or `NULL` for collapsed/union).
#' @param score_change_threshold Absolute score change used to mark `shifted`.
#' @param ... Additional arguments reserved for future extension.
#' @return A data.frame with at least `Feature`, `score`, `abs_score`, `sign`,
#'   `status`.
#' @noRd
extract_rewiring_features <- function(
    x,
    pathway_id,
    selected_layer,
    reference_layer = NULL,
    score_change_threshold = 0.1,
    ...
) {
    if (!inherits(x, "mnseaResult")) {
        stop("extract_rewiring_features() currently requires a mnseaResult.")
    }
    if (is.null(selected_layer)) {
        stop("selected_layer must be supplied.")
    }
    if (is.null(reference_layer)) {
        stop("reference_layer must be supplied; use the collapsed scores intentionally.")
    }

    selected_df <- fortify_mnsea_contribution(
        x,
        level = "feature",
        pathway_id = pathway_id,
        layer = selected_layer
    )
    reference_df <- fortify_mnsea_contribution(
        x,
        level = "feature",
        pathway_id = pathway_id,
        layer = reference_layer
    )

    selected_features <- setNames(selected_df$abs_score, selected_df$Feature)
    reference_features <- setNames(reference_df$abs_score, reference_df$Feature)

    all_features <- unique(c(names(selected_features), names(reference_features)))
    if (length(all_features) == 0) {
        return(data.frame(
            Feature = character(0),
            score = numeric(0),
            abs_score = numeric(0),
            sign = character(0),
            status = character(0),
            stringsAsFactors = FALSE
        ))
    }

    rows <- lapply(all_features, function(feature) {
        in_sel <- feature %in% names(selected_features)
        in_ref <- feature %in% names(reference_features)
        score <- if (in_sel) {
            selected_df$score[match(feature, selected_df$Feature)]
        } else {
            reference_df$score[match(feature, reference_df$Feature)]
        }
        abs_score <- if (in_sel) selected_features[[feature]] else reference_features[[feature]]
        sign <- if (score < 0) "suppressed" else if (score > 0) "activated" else "neutral"

        status <- if (in_sel && in_ref) {
            if (abs(abs_score - reference_features[[feature]]) > score_change_threshold) {
                "shifted"
            } else {
                "shared"
            }
        } else if (in_sel) {
            "gained"
        } else {
            "lost"
        }

        data.frame(
            Feature = feature,
            score = as.numeric(score),
            abs_score = as.numeric(abs_score),
            sign = sign,
            status = status,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

#' Classify mechanism state from enrichment shift and rewiring score
#'
#' @param nes_shift Numeric vector of NES differences.
#' @param rewiring_score Numeric vector of rewiring scores.
#' @param thresholds Named list with `rewire` and `nes` numeric thresholds.
#' @return Character vector of mechanism classes.
#' @noRd
classify_mechanism_state <- function(
    nes_shift,
    rewiring_score,
    thresholds = list(rewire = 0.5, nes = 0.2)
) {
    if (length(nes_shift) != length(rewiring_score)) {
        stop("nes_shift and rewiring_score must have the same length.")
    }
    vapply(seq_along(nes_shift), function(i) {
        d <- abs(nes_shift[i])
        r <- rewiring_score[i]
        if (is.na(r) || is.na(d)) {
            return(NA_character_)
        }
        if (r >= thresholds$rewire && d >= thresholds$nes) {
            "rewired"
        } else if (r >= thresholds$rewire) {
            "rewired"
        } else if (d >= thresholds$nes) {
            "context_specific"
        } else {
            "conserved"
        }
    }, character(1))
}

#' Summarize network mechanism information for a term-level table
#'
#' @param x A `nseaResult` or `mnseaResult`.
#' @param reference Optional reference result (another `nseaResult`,
#'   `mnseaResult`, or a results data.frame) used to compute `delta_NES`.
#' @param reference_layer Optional reference layer for mnsea comparisons.
#' @param selected_layer Optional layer to compare (defaults to collapsed).
#' @param thresholds Optional named list with `rewire` and `nes` numeric
#'   thresholds used by classify_mechanism_state().
#' @param ... Additional arguments passed to compute_rewiring_score().
#' @return A term-level data.frame with `reference_NES` and `delta_NES`.
#' @noRd
summarize_nsea_mechanism <- function(
    x,
    reference = NULL,
    reference_layer = NULL,
    selected_layer = NULL,
    thresholds = list(rewire = 0.5, nes = 0.2),
    ...
) {
    result_df <- .result_data(x)
    if (nrow(result_df) == 0) {
        return(result_df)
    }

    rew <- compute_rewiring_score(
        x,
        reference = reference,
        reference_layer = reference_layer,
        selected_layer = selected_layer,
        ...
    )

    ids <- as.character(result_df$ID)

    reference_NES <- rep(NA_real_, length(ids))
    delta_NES <- rep(NA_real_, length(ids))
    if (!is.null(reference)) {
        ref_df <- .result_data(reference)
        if (!"NES" %in% colnames(ref_df)) {
            stop("The reference result must contain an `NES` column.")
        }
        ref_NES <- setNames(as.numeric(ref_df$NES), as.character(ref_df$ID))
        reference_NES <- unname(ref_NES[ids])
        delta_NES <- as.numeric(result_df$NES) - reference_NES
    }

    if (nrow(rew) > 0 && all(ids %in% rew$ID)) {
        rew_match <- rew[match(ids, rew$ID), , drop = FALSE]
    } else {
        rew_match <- rew[0, , drop = FALSE]
    }

    leading_edge_size <- vapply(ids, function(id) {
        core <- result_df$core_enrichment[match(id, result_df$ID)]
        if (is.na(core) || !nzchar(core)) return(NA_integer_)
        length(strsplit(core, "/", fixed = TRUE)[[1]])
    }, integer(1))

    mechanism_class <- if (!is.null(reference) && nrow(rew_match) > 0) {
        classify_mechanism_state(
            nes_shift = delta_NES,
            rewiring_score = rew_match$rewiring_score,
            thresholds = thresholds
        )
    } else if (nrow(rew_match) > 0 && !all(is.na(rew_match$rewiring_score))) {
        # No reference NES is available; classify only on rewiring score.
        classify_mechanism_state(
            nes_shift = rep(0, nrow(rew_match)),
            rewiring_score = rew_match$rewiring_score,
            thresholds = thresholds
        )
    } else {
        rep(NA_character_, length(ids))
    }

    data.frame(
        ID = ids,
        Description = as.character(result_df$Description),
        NES = as.numeric(result_df$NES),
        reference_NES = reference_NES,
        delta_NES = delta_NES,
        p.adjust = as.numeric(result_df$p.adjust),
        leading_edge_size = leading_edge_size,
        leading_edge_overlap = if (nrow(rew_match) > 0) rew_match$leading_edge_overlap else NA_real_,
        rewiring_score = if (nrow(rew_match) > 0) rew_match$rewiring_score else NA_real_,
        centrality_shift = if (nrow(rew_match) > 0) rew_match$centrality_shift else NA_real_,
        mechanism_class = mechanism_class,
        stringsAsFactors = FALSE
    )
}

.mnsea_feature_sets <- function(x, layer = NULL) {
    ids <- as.character(x@result$ID)
    lapply(ids, function(id) {
        df <- fortify_mnsea_contribution(
            x,
            level = "feature",
            pathway_id = id,
            layer = layer
        )
        unique(as.character(df$Feature))
    }) |> stats::setNames(ids)
}

.result_feature_sets <- function(x, layer = NULL) {
    if (inherits(x, "mnseaResult")) {
        return(.mnsea_feature_sets(x, layer = layer))
    }
    if (inherits(x, "nseaResult")) {
        ids <- as.character(x@result$ID)
        gene_sets <- x@geneSets
        if (!is.list(gene_sets)) {
            return(stats::setNames(vector("list", length(ids)), ids))
        }
        return(lapply(ids, function(id) {
            set <- gene_sets[[id]]
            if (is.null(set)) character(0) else unique(as.character(set))
        }) |> stats::setNames(ids))
    }
    stop("x must be a nseaResult or mnseaResult.")
}
