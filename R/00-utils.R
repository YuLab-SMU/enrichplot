#' Internal package helpers loaded early
#'
#' @noRd
require_suggested <- function(package, feature) {
    for (pkg in as.character(package)) {
        rlang::check_installed(pkg, feature)
    }
    invisible(TRUE)
}

normalize_measure_var <- function(value, include_percentage = FALSE) {
    if (is.null(value) || length(value) != 1 || !is.character(value)) {
        return(value)
    }

    aliases <- c(
        geneRatio = "GeneRatio",
        GeneRatio = "GeneRatio",
        count = "Count",
        Count = "Count"
    )

    if (include_percentage) {
        aliases <- c(
            aliases,
            percentage = "Percentage",
            Percentage = "Percentage",
            rowPercentage = "Percentage"
        )
    }

    if (value %in% names(aliases)) {
        return(unname(aliases[[value]]))
    }

    value
}

compute_fold_enrichment <- function(df) {
    if (!all(c("GeneRatio", "BgRatio") %in% colnames(df))) {
        yulab.utils::yulab_abort(
            "`FoldEnrichment` is unavailable and cannot be derived without both `GeneRatio` and `BgRatio`."
        )
    }

    gene_ratio <- df$GeneRatio
    if (inherits(gene_ratio, "character")) {
        gene_ratio <- yulab.utils::parse_ratio(gene_ratio)
    }

    bg_ratio <- df$BgRatio
    if (inherits(bg_ratio, "character")) {
        bg_ratio <- yulab.utils::parse_ratio(bg_ratio)
    }

    as.numeric(gene_ratio) / as.numeric(bg_ratio)
}
