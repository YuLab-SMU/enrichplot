# Regression coverage for the tutorial-facing visualization functions.
#
# Guards the dispatch path reported by a reviewer (Comment 1.1, 2026-09
# revision round): under ggplot2 4.0.x with S7 < 0.2.2, `ggplot(df, aes(...))
# + theme_dose(...)` failed with `Incompatible methods ("Ops.S7_object",
# "+.gg") for "+"` and `non-numeric argument to binary operator`. Every
# plotting function below funnels through the same `+` dispatch, so a
# dispatch regression resurfaces in all of them.

make_rich_enrich_result <- function() {
    ids <- sprintf("T%d", 1:8)
    descriptions <- c(
        "cell cycle phase transition",
        "cell cycle DNA replication",
        "mitotic cell cycle process",
        "immune response regulation",
        "innate immune response",
        "cytokine signaling pathway",
        "lipid metabolic process",
        "fatty acid beta-oxidation"
    )
    gene_sets <- list(
        T1 = c("g1", "g2", "g3", "g4"),
        T2 = c("g1", "g2", "g3", "g5"),
        T3 = c("g1", "g2", "g4", "g6"),
        T4 = c("g5", "g6", "g7"),
        T5 = c("g5", "g7", "g8"),
        T6 = c("g6", "g7", "g9"),
        T7 = c("g8", "g9", "g10"),
        T8 = c("g8", "g10", "g1")
    )
    counts <- lengths(gene_sets)[ids]
    result <- data.frame(
        ID = ids,
        Description = descriptions,
        GeneRatio = paste0(counts, "/10"),
        BgRatio = paste0(counts + 20, "/500"),
        RichFactor = counts / 30,
        FoldEnrichment = counts / 30 * 2,
        zScore = seq(4, 1.1, length.out = 8),
        pvalue = c(1e-6, 3e-6, 1e-5, 1e-4, 2e-4, 5e-4, 1e-3, 4e-3),
        p.adjust = c(1e-5, 2e-5, 5e-5, 4e-4, 6e-4, 1.2e-3, 2e-3, 6e-3),
        qvalue = c(1e-5, 2e-5, 5e-5, 4e-4, 6e-4, 1.2e-3, 2e-3, 6e-3),
        geneID = vapply(gene_sets[ids], paste, character(1), collapse = "/"),
        Count = as.integer(counts),
        stringsAsFactors = FALSE
    )
    rownames(result) <- ids

    methods::new(
        "enrichResult",
        result = result,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2,
        organism = "mock",
        ontology = "mock",
        gene = paste0("g", 1:10),
        keytype = "UNKNOWN",
        universe = paste0("g", 1:500),
        gene2Symbol = character(),
        geneSets = gene_sets,
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

make_rich_gsea_result <- function() {
    ids <- sprintf("T%d", 1:4)
    descriptions <- c(
        "cell cycle phase transition",
        "immune response regulation",
        "lipid metabolic process",
        "fatty acid beta-oxidation"
    )
    gene_list <- c(
        g1 = 3, g2 = 2.5, g3 = 2, g4 = 1.5, g5 = 1.2,
        g6 = 0.8, g7 = 0.2, g8 = -0.5, g9 = -1, g10 = -2
    )
    gene_sets <- list(
        T1 = c("g1", "g2", "g3", "g4"),
        T2 = c("g4", "g5", "g6", "g7"),
        T3 = c("g6", "g7", "g8", "g9"),
        T4 = c("g8", "g9", "g10")
    )
    result <- data.frame(
        ID = ids,
        Description = descriptions,
        setSize = c(4L, 4L, 4L, 3L),
        enrichmentScore = c(0.85, 0.4, -0.6, -0.9),
        NES = c(1.8, 0.9, -1.2, -1.7),
        pvalue = c(0.01, 0.2, 0.1, 0.02),
        p.adjust = c(0.04, 0.5, 0.3, 0.08),
        qvalue = c(0.04, 0.5, 0.3, 0.08),
        rank = c(1L, 10L, 12L, 14L),
        leading_edge = c(
            "tags=100%, list=40%, signal=60%",
            "tags=50%, list=70%, signal=35%",
            "tags=75%, list=80%, signal=60%",
            "tags=100%, list=90%, signal=90%"
        ),
        core_enrichment = c(
            "g1/g2/g3/g4", "g4/g5", "g8/g9/g10", "g8/g9/g10"
        ),
        stringsAsFactors = FALSE
    )
    rownames(result) <- ids

    methods::new(
        "gseaResult",
        result = result,
        organism = "mock",
        setType = "mock",
        geneSets = gene_sets,
        geneList = gene_list,
        keytype = "UNKNOWN",
        permScores = matrix(runif(200), nrow = 50, ncol = 4),
        params = list(exponent = 1, nPerm = 50, pvalueCutoff = 0.25),
        gene2Symbol = character(),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

expect_ggplot <- function(p) {
    expect_s3_class(p, "ggplot")
    # force data-pipeline evaluation to catch missing-column / bad-aes errors
    if (inherits(p, "ggplot")) {
        expect_error(ggplot2::ggplot_build(p), NA)
    }
}

test_that("barplot.enrichResult runs the tutorial 13.1 call", {
    p <- barplot(make_rich_enrich_result(), showCategory = 5)
    expect_ggplot(p)
    expect_true(length(p$layers) > 0)
})

## ---------------------------------------------------------------------------
## The `+` dispatch canary: the exact expression shape that broke for the
## reviewer (ggplot2 4.0.x + S7 0.2.1). Keep independent of any mock object.
## ---------------------------------------------------------------------------

test_that("ggplot() + theme_dose() + geom dispatch resolves under ggplot2 4.x", {
    df <- data.frame(x = 1:3, y = c("a", "b", "c"), stringsAsFactors = FALSE)
    p <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]])) +
        theme_dose(12) +
        set_enrichplot_color(type = "fill", name = "p.adjust")
    p <- p + geom_col()
    expect_ggplot(p)
    expect_true(length(p$layers) > 0)
})

## ---------------------------------------------------------------------------
## barplot() — tutorial Section 13.1 entry point
## ---------------------------------------------------------------------------

test_that("barplot.enrichResult runs the tutorial 13.1 call", {
    p <- barplot(make_rich_enrich_result(), showCategory = 5)
    expect_ggplot(p)
    expect_true(length(p$layers) > 0)
})

test_that("barplot.enrichResult works when colorBy is absent from result", {
    x <- make_rich_enrich_result()
    x@result$zScore <- NULL
    expect_ggplot(barplot(x, showCategory = 5))
})

test_that("barplot methods forward width to geom_col", {
    p_enrich <- barplot(make_rich_enrich_result(), showCategory = 5, width = 0.2)
    p_compare <- barplot(mock_comparecluster_result(), showCategory = 2, width = 0.3)

    expect_equal(unique(ggplot2::ggplot_build(p_enrich)$data[[1]]$width), 0.2)
    expect_equal(unique(ggplot2::ggplot_build(p_compare)$data[[1]]$width), 0.3)
})

test_that("dotplot keeps hollow size legends after cowplot composition", {
    skip_if_not_installed("cowplot")

    p <- dotplot(make_rich_enrich_result(), showCategory = 5)
    combo <- cowplot::plot_grid(p, p, ncol = 1)
    grob <- grid::grid.force(ggplot2::ggplotGrob(combo))
    point_grobs <- grep("GRID.points", grid::grid.ls(grob, print = FALSE)$name, value = TRUE)
    legend_points <- lapply(point_grobs, function(name) grid::getGrob(grob, name, grep = TRUE))

    expect_true(length(legend_points) >= 6)
    expect_true(all(vapply(legend_points, `[[`, numeric(1), "pch") == enrichplot_point_shape))
    expect_true(all(vapply(legend_points, function(x) all(is.na(x$gp$fill)), logical(1))))
})

test_that("barplot.gseaResult runs", {
    expect_ggplot(barplot(make_rich_gsea_result(), showCategory = 3))
})

test_that("barplot.compareClusterResult runs", {
    expect_ggplot(barplot(mock_comparecluster_result(), showCategory = 2))
})

test_that("barplot.compareClusterResult runs for all documented by values", {
    # default by="geneRatio" used to crash with "object 'p' not found"
    # (no bar branch for it in plotting.clusterProfile); "rowPercentage" was
    # silently mis-mapped. Both go through the fortify-produced columns now.
    x <- mock_comparecluster_result()
    expect_ggplot(barplot(x, showCategory = 2, by = "geneRatio"))
    expect_ggplot(barplot(x, showCategory = 2, by = "count"))
    expect_ggplot(barplot(x, showCategory = 2, by = "rowPercentage"))
})

## ---------------------------------------------------------------------------
## dotplot()
## ---------------------------------------------------------------------------

test_that("dotplot methods run for enrichResult, gseaResult and compareClusterResult", {
    expect_ggplot(dotplot(make_rich_enrich_result(), showCategory = 5))
    expect_ggplot(dotplot(make_rich_gsea_result(), showCategory = 3))
    expect_ggplot(dotplot(mock_comparecluster_result(), showCategory = 2))
})

## ---------------------------------------------------------------------------
## Gene-concept network, heatmap, upset
## ---------------------------------------------------------------------------

test_that("cnetplot methods run", {
    expect_ggplot(cnetplot(make_rich_enrich_result(), showCategory = 3))
    expect_ggplot(cnetplot(make_rich_gsea_result(), showCategory = 2))
    expect_ggplot(cnetplot(mock_comparecluster_result(), showCategory = 2))
})

test_that("heatplot.enrichResult runs with and without foldChange", {
    x <- make_rich_enrich_result()
    expect_ggplot(heatplot(x, showCategory = 4))
    expect_ggplot(heatplot(x, showCategory = 4, foldChange = mock_foldchange()))
})

test_that("upsetplot methods run", {
    expect_ggplot(upsetplot(make_rich_enrich_result(), n = 4))
    expect_ggplot(upsetplot(make_rich_gsea_result(), n = 3))
})

test_that("upsetplot.gseaResult boxplots do not duplicate outliers", {
    p <- upsetplot(make_rich_gsea_result(), n = 3, type = "boxplot")

    expect_identical(p$layers[[1]]$geom_params$outlier_gp$shape, NA)
})

test_that("upsetplot.gseaResult remaps fold changes for readable objects", {
    x <- make_rich_gsea_result()
    x@readable <- TRUE
    x@gene2Symbol <- setNames(paste0("SYM", 1:10), paste0("g", 1:10))
    x@geneSets <- lapply(x@geneSets, function(gs) unname(x@gene2Symbol[gs]))
    x@result$core_enrichment <- vapply(
        x@geneSets[x@result$ID],
        paste,
        character(1),
        collapse = "/"
    )

    p <- upsetplot(x, n = 3)

    expect_false(anyNA(p$data$foldChange))
    expect_true(all(grepl("^SYM", p$data$gene)))
})

## ---------------------------------------------------------------------------
## Semantic-similarity based plots (emapplot / ssplot / treeplot)
## ---------------------------------------------------------------------------

test_that("pairwise_termsim fills the similarity matrix offline via JC", {
    x <- pairwise_termsim(make_rich_enrich_result(), method = "JC")
    expect_equal(nrow(x@termsim), 8)
    expect_equal(ncol(x@termsim), 8)
    expect_true(all(diag(as.matrix(x@termsim)) == 1))
})

test_that("emapplot / ssplot / treeplot run on JC similarity", {
    x <- pairwise_termsim(make_rich_enrich_result(), method = "JC")
    expect_ggplot(emapplot(x, showCategory = 8, nCluster = 2))
    expect_ggplot(ssplot(x, showCategory = 8))
    expect_ggplot(treeplot(x, showCategory = 8, nCluster = 2))
})

test_that("treeplot compareCluster heatMap panels run", {
    x <- pairwise_termsim(mock_comparecluster_result(), method = "JC")

    expect_ggplot(treeplot(x, showCategory = 2, cluster_panel = "heatMap"))
})

test_that("treeplot compareCluster dotplot panels run", {
    skip_if_not_installed("ggtreeExtra")

    x <- pairwise_termsim(mock_comparecluster_result(), method = "JC")

    expect_ggplot(treeplot(x, showCategory = 2, cluster_panel = "dotplot"))
})

test_that("emapplot and ssplot honor group_legend for grouped layouts", {
    x <- pairwise_termsim(make_rich_enrich_result(), method = "JC")

    p_emap_no_legend <- emapplot(
        x,
        showCategory = 8,
        group = TRUE,
        group_legend = FALSE
    )
    p_emap_with_legend <- emapplot(
        x,
        showCategory = 8,
        group = TRUE,
        group_legend = TRUE
    )
    p_ss <- ssplot(x, showCategory = 8, group_legend = FALSE)

    expect_ggplot(p_emap_no_legend)
    expect_ggplot(p_emap_with_legend)
    expect_ggplot(p_ss)
    expect_false(isTRUE(p_emap_no_legend$layers[[3]]$show.legend[["fill"]]))
    expect_true(isTRUE(p_emap_with_legend$layers[[3]]$show.legend[["fill"]]))
})

test_that("emapplot and ssplot accept label-keyed similarity matrices from non-JC methods", {
    # get_similarity_matrix() keys every method's matrix by term labels; a
    # regression re-mapped those keys as if they were IDs, producing NA edges
    # ("edge data frame contains NAs") for Wang and other semantic measures.
    x <- make_rich_enrich_result()
    labels <- x@result$Description
    sim <- matrix(0.5, nrow = 8, ncol = 8, dimnames = list(labels, labels))
    diag(sim) <- 1
    x@termsim <- sim
    x@method <- "Wang"
    expect_ggplot(emapplot(x, showCategory = 8))
    expect_ggplot(ssplot(x, showCategory = 8))
})

test_that("emapplot runs with real Wang similarity on GO terms", {
    skip_if_not_installed("GOSemSim")
    skip_if_not_installed("org.Hs.eg.db")
    go_ids <- c(
        "GO:0007049", "GO:0006260", "GO:0006952",
        "GO:0006629", "GO:0006631"
    )
    x <- make_rich_enrich_result()
    x@result <- x@result[1:5, ]
    x@result$ID <- go_ids
    x@result$Description <- c(
        "mitotic cell cycle",
        "DNA replication",
        "defense response",
        "lipid metabolic process",
        "fatty acid metabolic process"
    )
    rownames(x@result) <- go_ids
    x@ontology <- "BP"
    x@termsim <- matrix(0, 0, 0)
    sem_data <- suppressWarnings(GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP"))
    x <- pairwise_termsim(x, method = "Wang", semData = sem_data)
    expect_ggplot(emapplot(x, showCategory = 5))
    expect_ggplot(ssplot(x, showCategory = 5))
})

## ---------------------------------------------------------------------------
## Term-level statistical plots
## ---------------------------------------------------------------------------

test_that("volplot and manhattanplot run on enrichResult", {
    x <- make_rich_enrich_result()
    expect_ggplot(volplot(x, showCategory = 8))
    expect_ggplot(manhattanplot(x, showCategory = 8))
})

## ---------------------------------------------------------------------------
## GSEA plots
## ---------------------------------------------------------------------------

test_that("gseaplot methods run", {
    x <- make_rich_gsea_result()
    expect_ggplot(gseaplot(x, geneSetID = 1, by = "runningScore"))
    expect_ggplot(gseaplot(x, geneSetID = 1, by = "preranked"))
    expect_error(gseaplot(x, geneSetID = 1, by = "all"), NA)
})

test_that("gseaplot2 and gsearank run", {
    x <- make_rich_gsea_result()
    expect_error(gseaplot2(x, geneSetID = 1), NA)
    expect_error(gseaplot2(x, geneSetID = c(1, 2), subplots = 1), NA)
    expect_ggplot(gsearank(x, geneSetID = 1))
})

test_that("gseaplot2 gglist objects work with cowplot grids", {
    skip_if_not_installed("cowplot")

    x <- make_rich_gsea_result()
    p <- gseaplot2(
        x,
        geneSetID = 1,
        subplots = c(1, 2),
        pvalue_table = TRUE
    )

    expect_s3_class(p, "gglist")
    expect_s3_class(cowplot::as_grob(p), "grob")
    expect_error(cowplot::plot_grid(p, p, ncol = 1), NA)
})

test_that("gseaplot2 hit bins follow ranked-list order", {
    x <- make_rich_gsea_result()
    p <- gseaplot2(x, geneSetID = 1)
    rects <- ggplot2::ggplot_build(p[[2]])$data[[2]]

    expect_equal(rects$xmin, c(1, 2, 3, 4))
    expect_equal(rects$xmax, c(2, 3, 4, 11))
})

test_that("hplot runs", {
    expect_error(hplot(make_rich_gsea_result(), geneSetID = 1), NA)
})

test_that("gseadist runs", {
    expect_ggplot(gseadist(make_rich_gsea_result(), IDs = c(1, 2)))
})

test_that("ridgeplot.gseaResult runs", {
    skip_if_not_installed("ggridges")
    expect_ggplot(
        ridgeplot(
            make_rich_gsea_result(),
            showCategory = c("T1", "T3", "T4"),
            core_enrichment = TRUE
        )
    )
})

## ---------------------------------------------------------------------------
## Composition helpers
## ---------------------------------------------------------------------------

test_that("plot_list combines tutorial plots", {
    x <- make_rich_enrich_result()
    p <- plot_list(barplot(x, showCategory = 3), dotplot(x, showCategory = 3),
        ncol = 2, labels = c("A", "B")
    )
    expect_error(print(p), NA)
})
