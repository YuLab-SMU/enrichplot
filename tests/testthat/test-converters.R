test_that("import_enrichr maps enrichr table to enrichResult", {
    db <- data.frame(
        Term = c("GO:0006915;;apoptotic process", "Reactome Pathway X"),
        Overlap = c("3/50", "2/40"),
        P.value = c(0.001, 0.01),
        Adjusted.P.value = c(0.01, 0.05),
        Genes = c("G1; G2; G3", "G4;G5"),
        stringsAsFactors = FALSE
    )
    gene <- paste0("G", 1:10)
    universe <- paste0("G", 1:1000)
    suppressWarnings(
        x <- import_enrichr(db, gene = gene, universe = universe)
    )
    expect_s4_class(x, "enrichResult")
    expect_equal(x@result$ID, c("GO:0006915", "Reactome Pathway X"))
    expect_equal(x@result$Description, c("apoptotic process", "Reactome Pathway X"))
    expect_equal(x@result$GeneRatio, c("3/10", "2/10"))
    expect_equal(x@result$BgRatio, c("50/1000", "40/1000"))
    expect_equal(x@result$geneID, c("G1/G2/G3", "G4/G5"))
    expect_equal(x@result$Count, c(3, 2))
    expect_equal(x@result$RichFactor, c(3 / 50, 2 / 40))

    ## list input picks by name
    xl <- import_enrichr(
        list(KEGG = db, GO = db), db = "GO", gene = gene, universe = universe
    )
    expect_equal(xl@result$Count, c(3, 2))
    expect_error(import_enrichr(list(KEGG = db), db = "GO_Bio"), "not found")
})

test_that("import_gprofiler2 maps gost result", {
    gost_res <- list(result = data.frame(
        query.number = 1,
        precision = 0.3,
        query_size = 10,
        term_size = 50,
        effective_domain_size = 1000,
        term_id = "GO:0000001",
        term_name = "term one",
        p_value = 0.001,
        adjusted_p_value = 0.01,
        intersection_size = 3L,
        intersections = I(list(c("G1", "G2", "G3"))),
        stringsAsFactors = FALSE
    ))
    x <- suppressWarnings(import_gprofiler2(gost_res))
    expect_s4_class(x, "enrichResult")
    expect_equal(x@result$ID, "GO:0000001")
    expect_equal(x@result$GeneRatio, "3/10")
    expect_equal(x@result$BgRatio, "50/1000")
    expect_equal(x@result$geneID, "G1/G2/G3")
    expect_equal(x@result$p.adjust, 0.01)
    expect_equal(x@result$Count, 3)

    ## falls back to 'gene' when query_size is absent
    gost_res$result$query_size <- NULL
    x2 <- suppressWarnings(
        import_gprofiler2(gost_res, gene = paste0("G", 1:10))
    )
    expect_equal(x2@result$GeneRatio, "3/10")
})

test_that("import_webgestalt maps WebGestaltR table", {
    wg <- data.frame(
        geneSet = "hsa04110",
        description = "cell cycle",
        size = 50L,
        overlap = 3L,
        rawPValue = 0.001,
        adjPValue = 0.01,
        userIds = "G1;G2;G3",
        stringsAsFactors = FALSE
    )
    x <- suppressWarnings(
        import_webgestalt(wg, gene = paste0("G", 1:10), universe = paste0("G", 1:1000))
    )
    expect_s4_class(x, "enrichResult")
    expect_equal(x@result$ID, "hsa04110")
    expect_equal(x@result$Description, "cell cycle")
    expect_equal(x@result$GeneRatio, "3/10")
    expect_equal(x@result$BgRatio, "50/1000")
    expect_equal(x@result$geneID, "G1/G2/G3")
    expect_equal(x@result$RichFactor, 3 / 50)
})

test_that("import_fgsea maps fgsea result to gseaResult", {
    stats <- c(G1 = 3, G2 = 2, G3 = 1, G4 = -1, G5 = -2, G6 = -3)
    pathways <- list(P1 = c("G1", "G2", "G3"), P2 = c("G4", "G5", "G6"))
    fgres <- data.frame(
        pathway = c("P1", "P2"),
        ES = c(0.6, -0.5),
        NES = c(1.5, -1.4),
        pval = c(0.01, 0.02),
        padj = c(0.02, 0.02),
        size = c(3L, 3L),
        stringsAsFactors = FALSE
    )
    fgres$leadingEdge <- list(c("G1", "G2"), c("G5", "G6"))
    x <- suppressWarnings(import_fgsea(fgres, stats = stats, geneSets = pathways))

    expect_s4_class(x, "gseaResult")
    expect_equal(x@geneList, stats)
    expect_equal(names(x@geneSets), c("P1", "P2"))
    expect_true(all(c("rank", "leading_edge", "core_enrichment", "setSize") %in%
        colnames(x@result)))
    expect_equal(x@result$core_enrichment[1], "G1/G2")
    expect_true(all(x@result$rank > 0))
    expect_equal(x@result$setSize, c(3, 3))

    ## without geneSets, rebuilt from leading edge with a warning
    expect_warning(
        x2 <- import_fgsea(fgres, stats = stats),
        "geneSets"
    )
    expect_s4_class(x2, "gseaResult")
})

test_that("importers give informative errors on missing columns", {
    expect_error(import_enrichr(data.frame(Term = "X")), "missing required column")
    expect_error(
        import_gprofiler2(data.frame(term_id = "X")),
        "missing required column"
    )
    expect_error(
        import_webgestalt(data.frame(geneSet = "P")),
        "missing required column"
    )
    expect_error(
        import_fgsea(data.frame(pathway = "P", ES = 0.5)),
        "missing required column"
    )
})

test_that("imported results work with enrichplot visualization", {
    db <- data.frame(
        Term = c("GO:0006915;;apoptotic process", "Reactome Pathway X"),
        Overlap = c("3/50", "2/40"),
        P.value = c(0.001, 0.01),
        Adjusted.P.value = c(0.01, 0.05),
        Genes = c("G1; G2; G3", "G4;G5"),
        stringsAsFactors = FALSE
    )
    x <- suppressWarnings(
        import_enrichr(db, gene = paste0("G", 1:10), universe = paste0("G", 1:1000))
    )
    p <- dotplot(x)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    p2 <- barplot(x)
    expect_s3_class(ggplot2::ggplot_build(p2), "ggplot_built")

    stats <- c(G1 = 3, G2 = 2, G3 = 1, G4 = -1, G5 = -2, G6 = -3)
    pathways <- list(P1 = c("G1", "G2", "G3"), P2 = c("G4", "G5", "G6"))
    fgres <- data.frame(
        pathway = c("P1", "P2"),
        ES = c(0.6, -0.5),
        NES = c(1.5, -1.4),
        pval = c(0.01, 0.02),
        padj = c(0.02, 0.02),
        size = c(3L, 3L),
        stringsAsFactors = FALSE
    )
    fgres$leadingEdge <- list(c("G1", "G2"), c("G5", "G6"))
    gx <- suppressWarnings(import_fgsea(fgres, stats = stats, geneSets = pathways))
    pg <- gseaplot(gx, geneSetID = "P1")
    plots <- if (inherits(pg, "gglist")) pg else list(pg)
    built <- lapply(plots, ggplot2::ggplot_build)
    expect_true(all(vapply(built, inherits, logical(1), "ggplot_built")))
})

test_that("external importers do not claim a clusterProfiler citation", {
    db <- data.frame(
        Term = "GO:0006915;;apoptotic process",
        Overlap = "1/10",
        P.value = 0.001,
        Adjusted.P.value = 0.01,
        Genes = "G1",
        stringsAsFactors = FALSE
    )
    ora <- suppressWarnings(import_enrichr(db, gene = "G1"))

    fgsea_result <- data.frame(
        pathway = "P1", ES = 0.5, pval = 0.01,
        stringsAsFactors = FALSE
    )
    gsea <- suppressWarnings(
        import_fgsea(
            fgsea_result,
            stats = c(G1 = 1),
            geneSets = list(P1 = "G1")
        )
    )

    ora_output <- capture.output(show(ora))
    gsea_output <- capture.output(show(gsea))
    expect_false(any(grepl("clusterProfiler", ora_output, fixed = TRUE)))
    expect_false(any(grepl("clusterProfiler", gsea_output, fixed = TRUE)))
})
