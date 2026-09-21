expect_ggplot_build_ok <- function(p) {
    expect_s3_class(p, "ggplot")
    expect_error(ggplot2::ggplot_build(p), NA)
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

test_that("pairwise_termsim uses stable labels when descriptions repeat", {
    x <- mock_enrich_result()

    y <- pairwise_termsim(x, method = "JC", showCategory = c("T1", "T2"))

    expect_equal(rownames(y@termsim), c("dup [T1]", "dup [T2]"))
    expect_equal(colnames(y@termsim), c("dup [T1]", "dup [T2]"))
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
