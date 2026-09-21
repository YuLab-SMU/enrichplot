expect_ggplot_build_ok <- function(p) {
    expect_s3_class(p, "ggplot")
    expect_error(ggplot2::ggplot_build(p), NA)
}

make_strict_cutoff_enrich_result <- function() {
    result <- data.frame(
        ID = c("T1", "T2"),
        Description = c("a", "b"),
        GeneRatio = c("1/2", "1/2"),
        BgRatio = c("1/10", "1/10"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.2, 0.3),
        qvalue = c(0.2, 0.3),
        geneID = c("g1/g2", "g2/g3"),
        Count = c(2L, 2L),
        stringsAsFactors = FALSE
    )
    rownames(result) <- result$ID

    methods::new(
        "enrichResult",
        result = result,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        organism = "mock",
        ontology = "mock",
        gene = c("g1", "g2", "g3"),
        keytype = "UNKNOWN",
        universe = character(),
        gene2Symbol = character(),
        geneSets = list(
            T1 = c("g1", "g2"),
            T2 = c("g2", "g3")
        ),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

make_short_core_gsea_result <- function() {
    result <- data.frame(
        ID = c("T1", "T2", "T3"),
        Description = c("dense core set", "short core set A", "short core set B"),
        setSize = c(3L, 2L, 2L),
        enrichmentScore = c(0.9, 0.7, 0.6),
        NES = c(1.8, 1.4, 1.2),
        pvalue = c(0.01, 0.02, 0.03),
        p.adjust = c(0.02, 0.03, 0.04),
        qvalue = c(0.02, 0.03, 0.04),
        rank = c(1L, 2L, 3L),
        leading_edge = c(
            "tags=100%, list=40%, signal=60%",
            "tags=100%, list=60%, signal=50%",
            "tags=100%, list=80%, signal=40%"
        ),
        core_enrichment = c("g1/g2/g3", "g4/g5", "g4/g5"),
        stringsAsFactors = FALSE
    )
    rownames(result) <- result$ID

    methods::new(
        "gseaResult",
        result = result,
        organism = "mock",
        setType = "mock",
        geneSets = list(
            T1 = c("g1", "g2", "g3"),
            T2 = c("g4", "g5"),
            T3 = c("g4", "g5")
        ),
        geneList = c(g1 = 3, g2 = 2, g3 = 1, g4 = -1, g5 = -2),
        keytype = "UNKNOWN",
        permScores = matrix(runif(30), nrow = 10, ncol = 3),
        params = list(exponent = 1, nPerm = 10, pvalueCutoff = 1),
        gene2Symbol = character(),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

test_that("duplicate term descriptions get stable display labels", {
    x <- mock_enrich_result()

    mapping <- enrichplot:::get_term_mapping(x)

    expect_equal(mapping$ID, c("T1", "T2"))
    expect_equal(mapping$label, c("dup [T1]", "dup [T2]"))
})

test_that("extract_geneSets keeps ID selection and carries display labels", {
    x <- mock_enrich_result()

    gene_sets <- enrichplot:::extract_geneSets(x, c("T2", "dup [T1]"))

    expect_equal(names(gene_sets), c("T2", "T1"))
    expect_equal(
        attr(gene_sets, "term_labels"),
        c(T2 = "dup [T2]", T1 = "dup [T1]")
    )
})

test_that("heatplot(showTop) fails early without foldChange", {
    x <- mock_enrich_result()

    expect_error(
        heatplot(x, showCategory = c("T1", "T2"), showTop = 1),
        "`showTop` requires `foldChange`."
    )
})

test_that("heatplot uses disambiguated labels for duplicate descriptions", {
    x <- mock_enrich_result()

    p <- heatplot(
        x,
        showCategory = c("T1", "T2"),
        showTop = 1,
        foldChange = mock_foldchange()
    )

    expect_s3_class(p, "ggplot")
    expect_setequal(unique(p$data$categoryID), c("dup [T1]", "dup [T2]"))
    expect_equal(
        anyDuplicated(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()),
        0L
    )
})

test_that("ridgeplot drops undersized core gene sets instead of drawing empty rows", {
    skip_if_not_installed("ggridges")

    x <- make_short_core_gsea_result()

    expect_warning(
        p <- ridgeplot(x, showCategory = 3, core_enrichment = TRUE),
        "Dropping 2 gene set\\(s\\) with fewer than 3 ranked values"
    )

    expect_ggplot_build_ok(p)
    expect_equal(unique(as.character(p$data$category)), "dense core set")
    expect_equal(
        ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels(),
        "dense core set"
    )
})

test_that("pairwise_termsim uses stable labels when descriptions repeat", {
    x <- mock_enrich_result()

    y <- pairwise_termsim(x, method = "JC", showCategory = c("T1", "T2"))

    expect_equal(rownames(y@termsim), c("dup [T1]", "dup [T2]"))
    expect_equal(colnames(y@termsim), c("dup [T1]", "dup [T2]"))
})

test_that("pairwise_termsim can use raw enrichResult rows beyond significance cutoffs", {
    x <- make_strict_cutoff_enrich_result()

    expect_equal(nrow(as.data.frame(x)), 0)

    y <- pairwise_termsim(x, method = "JC", showCategory = 2)

    expect_equal(dim(y@termsim), c(2L, 2L))
    expect_equal(rownames(y@termsim), c("a", "b"))
    expect_ggplot_build_ok(emapplot(y, showCategory = 2))
})

test_that("get_enrichplot_color expands two custom colors to three safely", {
    old <- getOption("enrichplot.colours")
    on.exit(options(enrichplot.colours = old), add = TRUE)

    options(enrichplot.colours = c("#111111", "#222222"))

    expect_equal(
        enrichplot:::get_enrichplot_color(3),
        c("#111111", "white", "#222222")
    )
})

test_that("cnetplot smoke test works for compareClusterResult", {
    x <- mock_comparecluster_result()

    p <- cnetplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
})

test_that("cnetplot keeps duplicate compareCluster descriptions as separate term nodes", {
    x <- mock_comparecluster_result()
    d <- as.data.frame(x@compareClusterResult, stringsAsFactors = FALSE)
    d$ID <- c("T1", "T2", "T3", "T4")
    d$Description <- c("dup", "other", "dup", "other")
    d$geneID <- c("1/2", "2/3", "4/5", "5/6")
    x@compareClusterResult <- d

    p <- cnetplot(x, showCategory = 4)

    expect_ggplot_build_ok(p)
    expect_setequal(
        unique(as.character(p$data$name[p$data$.isCategory])),
        c("dup [T1]", "other [T2]", "dup [T3]", "other [T4]")
    )
})

test_that("emapplot smoke test works for compareClusterResult", {
    x <- mock_comparecluster_result()
    x <- pairwise_termsim(x, method = "JC", showCategory = 2)

    p <- emapplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
})

test_that("compareCluster pie plots tolerate duplicated cluster-term rows", {
    x <- mock_comparecluster_result()
    d <- x@compareClusterResult[c(1, 3, 2), , drop = FALSE]
    d$Cluster <- factor(c("A", "A", "A"), levels = c("A", "B"))
    x@compareClusterResult <- d

    expect_ggplot_build_ok(cnetplot(x, showCategory = 3))

    x <- pairwise_termsim(x, method = "JC", showCategory = 3)
    expect_ggplot_build_ok(emapplot(x, showCategory = 3))
})

test_that("manhattanplot normalizes lowercase size aliases", {
    expect_ggplot_build_ok(
        manhattanplot(mock_enrich_result(), showCategory = 2, size = "count")
    )
    expect_ggplot_build_ok(
        manhattanplot(
            mock_comparecluster_result(),
            showCategory = 2,
            size = "geneRatio"
        )
    )
})

test_that("dotplot supports Percentage as a size measure", {
    p <- dotplot(
        mock_enrich_result(),
        x = "geneRatio",
        size = "Percentage",
        showCategory = 2
    )

    expect_ggplot_build_ok(p)
    expect_true("Percentage" %in% colnames(p$data))
    expect_equal(p$data$Percentage, p$data$GeneRatio * 100)
})

test_that("dotplot2 derives FoldEnrichment from ratio columns", {
    p <- dotplot2(
        mock_comparecluster_result(),
        vars = c("A", "B"),
        label = c(A = "Control", B = "Treated")
    )

    expect_ggplot_build_ok(p)
})

test_that("dotplot2 fails clearly when x cannot be derived", {
    x <- mock_comparecluster_result()
    x@compareClusterResult$BgRatio <- NULL

    expect_error(
        dotplot2(x, vars = c("A", "B")),
        "`FoldEnrichment` is unavailable"
    )
})

test_that("dotplot applies numeric showCategory after orderBy sorting", {
    x <- methods::new(
        "enrichResult",
        result = structure(
            data.frame(
                ID = paste0("T", 1:4),
                Description = paste0("term", 1:4),
                GeneRatio = c("1/10", "2/10", "3/10", "4/10"),
                BgRatio = rep("1/100", 4),
                pvalue = c(0.05, 0.001, 0.02, 0.03),
                p.adjust = c(0.05, 0.001, 0.02, 0.03),
                qvalue = c(0.05, 0.001, 0.02, 0.03),
                geneID = c("g1", "g1/g2", "g1/g2/g3", "g1/g2/g3/g4"),
                Count = c(1L, 50L, 10L, 20L),
                stringsAsFactors = FALSE
            ),
            row.names = paste0("T", 1:4)
        ),
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        organism = "mock",
        ontology = "mock",
        gene = c("g1", "g2", "g3", "g4"),
        keytype = "UNKNOWN",
        universe = character(),
        gene2Symbol = character(),
        geneSets = list(
            T1 = "g1",
            T2 = c("g1", "g2"),
            T3 = c("g1", "g2", "g3"),
            T4 = c("g1", "g2", "g3", "g4")
        ),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )

    ids2 <- ggplot2::ggplot_build(
        dotplot(x, showCategory = 2, orderBy = "Count")
    )$plot$data$ID
    ids4 <- ggplot2::ggplot_build(
        dotplot(x, showCategory = 4, orderBy = "Count")
    )$plot$data$ID

    expect_identical(ids2, c("T2", "T4"))
    expect_identical(ids4[seq_along(ids2)], ids2)
})
