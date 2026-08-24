test_that("fortify_mnsea_contribution standardizes pathway tables", {
    x <- mock_mnsea_result()

    df <- fortify_mnsea_contribution(x, level = "pathway")

    expect_true(all(
        c("ID", "Description", "layer", "contribution", "share", "n_feature") %in%
            colnames(df)
    ))
    expect_true(nrow(df) > 0)
})

test_that("fortify_mnsea_subnetwork standardizes nodes and edges", {
    x <- mock_mnsea_result()

    subnet <- fortify_mnsea_subnetwork(x, pathway_id = "T1")

    expect_true(all(
        c("node", "label", "node_type", "layer", "abs_score", "sign") %in%
            colnames(subnet$nodes)
    ))
    expect_true(all(
        c("from", "to", "edge_type", "abs_weight") %in%
            colnames(subnet$edges)
    ))
    expect_true("pathway" %in% subnet$nodes$node_type)
    expect_true(any(subnet$edges$edge_type %in% c("membership", "intra", "coupling")))
})

test_that("mnsea helpers use a stable default pathway when pathway_id is NULL", {
    x <- mock_mnsea_result()

    feature_df <- fortify_mnsea_contribution(x, level = "feature")
    subnet <- fortify_mnsea_subnetwork(x)

    expect_true(all(feature_df$ID == "T1"))
    expect_true(all(subnet$nodes$pathway_id == "T1"))
    expect_true(all(subnet$pathway$ID == "T1"))
})

test_that("mnsea helpers support stable layer filtering", {
    x <- mock_mnsea_result()

    pathway_df <- fortify_mnsea_contribution(x, level = "pathway", layer = "rna")
    feature_df <- fortify_mnsea_contribution(x, level = "feature", pathway_id = "T1", layer = "protein")
    subnet <- fortify_mnsea_subnetwork(x, pathway_id = "T1", layer = "rna")

    expect_true(all(as.character(pathway_df$layer) == "rna"))
    expect_true(all(as.character(feature_df$layer) == "protein"))
    expect_true(all(as.character(subnet$nodes$layer[subnet$nodes$node_type == "feature"]) == "rna"))
})

test_that("dotplot works for mnseaResult contribution view", {
    x <- mock_mnsea_result()

    p <- dotplot(
        x,
        x = "share",
        color = "contribution",
        showCategory = 2
    )

    expect_s3_class(p, "ggplot")
    expect_true(all(c("layer", "share", "contribution") %in% colnames(p$data)))
})

test_that("dotplot mnseaResult uses layer filter and readable legend names", {
    x <- mock_mnsea_result()

    p <- dotplot(
        x,
        x = "share",
        color = "contribution",
        showCategory = 2,
        layer = "rna"
    )

    expect_true(all(as.character(p$data$layer) == "rna"))
    expect_equal(p$scales$get_scales("fill")$name, "Contribution")
    expect_equal(p$scales$get_scales("size")$name, "Feature count")
})

test_that("heatplot works for mnseaResult pathway contribution view", {
    x <- mock_mnsea_result()

    p <- heatplot(x, showCategory = 2, value = "contribution")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("Description", "layer", "fill_value") %in% colnames(p$data)))
    expect_setequal(unique(as.character(p$data$Description)), c("Pathway 1", "Pathway 2"))
})

test_that("heatplot mnseaResult uses layer filter and readable fill labels", {
    x <- mock_mnsea_result()

    p <- heatplot(x, showCategory = 2, layer = "protein", value = "contribution")

    expect_true(all(as.character(p$data$layer) == "protein"))
    expect_equal(p$scales$get_scales("fill")$name, "Contribution")
})

test_that("emapplot works for mnseaResult without precomputed termsim", {
    x <- mock_mnsea_result()

    p <- emapplot(x, showCategory = 2, min_edge = 0)

    expect_s3_class(p, "ggplot")
    expect_true(all(c("label", "p.adjust") %in% colnames(p$data)))
    expect_setequal(unique(as.character(p$data$label)), c("Pathway 1", "Pathway 2"))
})

test_that("emapplot mnseaResult supports layer filter and readable legends", {
    x <- mock_mnsea_result()

    p <- emapplot(x, showCategory = 2, layer = "rna", color = "contribution", min_edge = 0)

    expect_s3_class(p, "ggplot")
    expect_equal(unname(p$data$contribution[match(c("Pathway 1", "Pathway 2"), p$data$label)]), c(0.7, 0.4))
    expect_equal(p$scales$get_scales("colour")$name, "Contribution")
    expect_equal(p$scales$get_scales("size")$name, "Feature count")
})

test_that("emapplot mnseaResult respects stable term selection labels", {
    x <- mock_mnsea_result()

    p <- emapplot(x, showCategory = c("T2", "Pathway 1"), min_edge = 0)

    expect_s3_class(p, "ggplot")
    expect_setequal(unique(as.character(p$data$label)), c("Pathway 1", "Pathway 2"))
})

test_that("heatplot works for mnseaResult pathway-specific feature view", {
    x <- mock_mnsea_result()

    p <- heatplot(x, pathway_id = "T1", showTop = 2, value = "score")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("Feature", "layer", "fill_value") %in% colnames(p$data)))
    expect_lte(length(unique(as.character(p$data$Feature))), 2)
    expect_true(all(p$data$ID == "T1"))
})

test_that("gseaplot works for mnseaResult with collapsed scores", {
    x <- mock_mnsea_result()

    p <- gseaplot(x, geneSetID = "T1", by = "runningScore")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("ID", "Description", "runningScore", "position", "layer") %in% colnames(p$data)))
    expect_equal(unique(as.character(p$data$layer)), "collapsed")
    expect_equal(unique(as.character(p$data$ID)), "T1")
    expect_equal(as.numeric(p$data$geneList), as.numeric(x@collapsed_scores))
})

test_that("gseaplot mnseaResult supports single-layer ranked scores", {
    x <- mock_mnsea_result()

    p <- gseaplot(x, geneSetID = "T1", by = "preranked", layer = "protein")

    expect_s3_class(p, "ggplot")
    expect_equal(unique(as.character(p$data$layer)), "protein")
    expect_equal(as.numeric(p$data$geneList), as.numeric(x@layer_scores$protein))
})

test_that("gseaplot mnseaResult resolves numeric geneSetID to stable pathway IDs", {
    x <- mock_mnsea_result()

    p_idx <- gseaplot(x, geneSetID = 1, by = "runningScore")
    p_id <- gseaplot(x, geneSetID = "T1", by = "runningScore")

    expect_equal(as.numeric(p_idx$data$runningScore), as.numeric(p_id$data$runningScore))
    expect_equal(as.character(p_idx$data$Description), as.character(p_id$data$Description))
})

test_that("gseaplot mnseaResult rejects multi-layer requests", {
    x <- mock_mnsea_result()

    expect_error(
        gseaplot(x, geneSetID = "T1", by = "runningScore", layer = c("rna", "protein")),
        "`layer` must be `NULL` or a single layer"
    )
})

test_that("ridgeplot works for mnseaResult with collapsed scores", {
    x <- mock_mnsea_result()

    p <- ridgeplot(x, showCategory = 2, fill = "NES")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("ID", "category", "value", "layer", "NES") %in% colnames(p$data)))
    expect_equal(unique(as.character(p$data$layer)), "collapsed")
    expect_setequal(unique(as.character(p$data$ID)), c("T1", "T2"))
})

test_that("ridgeplot mnseaResult supports single-layer ranked scores", {
    x <- mock_mnsea_result()

    p <- ridgeplot(x, showCategory = "T1", fill = "NES", layer = "protein")

    expect_s3_class(p, "ggplot")
    expect_equal(unique(as.character(p$data$layer)), "protein")
    expect_equal(unique(as.character(p$data$ID)), "T1")
    expect_setequal(unique(as.numeric(p$data$value)), c(-0.4, 1.0))
})

test_that("ridgeplot mnseaResult respects core enrichment filter", {
    x <- mock_mnsea_result()

    p_core <- ridgeplot(
        x,
        showCategory = "T2",
        fill = "NES",
        core_enrichment = TRUE,
        layer = "rna"
    )
    p_all <- ridgeplot(
        x,
        showCategory = "T2",
        fill = "NES",
        core_enrichment = FALSE,
        layer = "rna"
    )

    expect_false(-0.5 %in% as.numeric(p_core$data$value))
    expect_true(-0.5 %in% as.numeric(p_all$data$value))
})

test_that("ridgeplot mnseaResult rejects multi-layer requests", {
    x <- mock_mnsea_result()

    expect_error(
        ridgeplot(x, showCategory = "T1", fill = "NES", layer = c("rna", "protein")),
        "`layer` must be `NULL` or a single layer"
    )
})

test_that("upsetplot works for mnseaResult with collapsed scores", {
    x <- mock_mnsea_result()
    df <- build_mnsea_upset_df(x, n = 2)
    p <- upsetplot(x, n = 2)

    expect_s3_class(p, "ggplot")
    expect_true(all(c("Feature", "score", "abs_score", "layer", "Description") %in% colnames(df)))
    expect_equal(unique(as.character(df$layer)), "collapsed")
    expect_setequal(unique(as.character(df$Feature)), c("g1", "g2", "g3"))
})

test_that("upsetplot mnseaResult supports single-layer feature overlap", {
    x <- mock_mnsea_result()
    df <- build_mnsea_upset_df(x, n = 2, layer = "protein", value = "abs_score")
    p <- upsetplot(x, n = 2, layer = "protein", value = "abs_score")

    expect_s3_class(p, "ggplot")
    expect_equal(unique(as.character(df$layer)), "protein")
    expect_setequal(unique(as.numeric(df$abs_score)), c(0.4, 0.9, 1.0))
})

test_that("upsetplot mnseaResult respects core enrichment filter", {
    x <- mock_mnsea_result()
    df_core <- build_mnsea_upset_df(x, n = 2, layer = "rna", core_enrichment = TRUE)
    df_all <- build_mnsea_upset_df(x, n = 2, layer = "rna", core_enrichment = FALSE)

    expect_false("g3" %in% as.character(df_core$Feature))
    expect_true("g3" %in% as.character(df_all$Feature))
})

test_that("ssplot works for mnseaResult with collapsed similarity space", {
    x <- mock_mnsea_result()

    p <- ssplot(x, showCategory = 2, node_label = "none")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("label", "size") %in% colnames(p$data)))
    expect_setequal(unique(as.character(p$data$label)), c("Pathway 1", "Pathway 2"))
})

test_that("ssplot mnseaResult supports single-layer similarity space", {
    x <- mock_mnsea_result()
    dr <- get_mnsea_drResult(x, showCategory = 2, layer = "rna", drfun = stats::cmdscale, dr.params = list(eig = TRUE))

    expect_true(all(c("Dimension1", "Dimension2") %in% colnames(dr$drdata)))
    expect_setequal(rownames(dr$drdata), c("Pathway 1", "Pathway 2"))
})

test_that("cnetplot uses the default mnsea pathway when pathway_id is NULL", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, node_label = "category")

    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$pathway_id == "T1"))
})

test_that("cnetplot mnseaResult supports layer filtering", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", layer = "protein", node_label = "none")

    feature_layers <- unique(as.character(p$data$layer[p$data$node_type == "feature"]))
    expect_equal(feature_layers, "protein")
})

test_that("heatplot uses the default mnsea pathway for feature view", {
    x <- mock_mnsea_result()

    p <- heatplot(x, showTop = 2, value = "score")

    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$ID == "T1"))
    expect_lte(length(unique(as.character(p$data$Feature))), 2)
})

test_that("heatplot rejects incompatible mnseaResult value semantics", {
    x <- mock_mnsea_result()

    expect_error(
        heatplot(x, pathway_id = "T1", value = "share"),
        "When `pathway_id` is provided, `value` must be `score` or `abs_score`."
    )
})

test_that("cnetplot works for mnseaResult subnetworks", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", include_couplings = TRUE, node_label = "category")

    expect_s3_class(p, "ggplot")
    expect_true("pathway" %in% p$data$node_type)
})

test_that("cnetplot supports share and exclusive mnsea labels", {
    x <- mock_mnsea_result()

    p_share <- cnetplot(x, pathway_id = "T1", node_label = "share")
    p_exclusive <- cnetplot(x, pathway_id = "T1", node_label = "exclusive")

    share_label_layers <- Filter(function(layer) inherits(layer$geom, "GeomTextRepel"), p_share$layers)
    exclusive_label_layers <- Filter(function(layer) inherits(layer$geom, "GeomTextRepel"), p_exclusive$layers)
    share_labels <- unique(unlist(lapply(share_label_layers, function(layer) layer$data$label)))
    exclusive_labels <- unique(unlist(lapply(exclusive_label_layers, function(layer) layer$data$label)))

    expect_true("g1" %in% share_labels || "g2" %in% share_labels)
    expect_true(length(exclusive_labels) == 0 || "g1" %in% exclusive_labels || "g2" %in% exclusive_labels)
})

test_that("cnetplot default mnsea labels keep pathway labels and deduplicate features", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "all")
    label_layers <- Filter(function(layer) inherits(layer$geom, "GeomTextRepel"), p$layers)
    label_data <- do.call(rbind, lapply(label_layers, function(layer) layer$data))

    expect_true("Pathway 1" %in% label_data$label)
    feature_labels <- label_data$label[label_data$node_type == "feature"]
    expect_equal(anyDuplicated(feature_labels), 0L)
})

test_that("cnetplot item labels use one representative label per feature", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "item")
    label_layers <- Filter(function(layer) inherits(layer$geom, "GeomTextRepel"), p$layers)
    item_labels <- unique(unlist(lapply(label_layers, function(layer) layer$data$label)))

    expect_equal(anyDuplicated(item_labels), 0L)
})

test_that("mnsea feature label selection prefers shared features when capped", {
    feature_nodes <- data.frame(
        Feature = c("g1", "g1", "g2", "g3"),
        plot_size = c(0.6, 0.5, 1.5, 1.2),
        membership_class = c("share", "share", "exclusive", "exclusive"),
        stringsAsFactors = FALSE
    )

    selected <- enrichplot:::select_mnsea_feature_labels(feature_nodes, max_labels = 2)

    expect_equal(selected$Feature, c("g1", "g2"))
})

test_that("mnsea label data splits pathway and feature labels for all mode", {
    pathway_nodes <- data.frame(
        label = "Pathway 1",
        node_type = "pathway",
        stringsAsFactors = FALSE
    )
    feature_nodes <- data.frame(
        label = c("g1", "g1", "g2"),
        Feature = c("g1", "g1", "g2"),
        plot_size = c(1.2, 0.8, 0.7),
        membership_class = c("share", "share", "exclusive"),
        node_type = "feature",
        stringsAsFactors = FALSE
    )

    label_data <- enrichplot:::select_mnsea_label_data(
        node_label = "all",
        node_data = feature_nodes[0, , drop = FALSE],
        pathway_nodes = pathway_nodes,
        feature_nodes = feature_nodes
    )

    expect_equal(label_data$pathway$label, "Pathway 1")
    expect_equal(label_data$feature$label, c("g1", "g2"))
})

test_that("cnetplot encodes mnsea edge semantics and honors size_edge", {
    x <- mock_mnsea_result()

    p_small <- cnetplot(x, pathway_id = "T1", include_couplings = TRUE, size_edge = 0.5)
    p_large <- cnetplot(x, pathway_id = "T1", include_couplings = TRUE, size_edge = 2)
    p_no_coupling <- cnetplot(x, pathway_id = "T1", include_couplings = FALSE, size_edge = 2)

    edge_small <- ggplot2::ggplot_build(p_small)$data[[1]]
    edge_large <- ggplot2::ggplot_build(p_large)$data[[1]]
    edge_no_coupling <- ggplot2::ggplot_build(p_no_coupling)$data[[1]]

    expect_gt(mean(edge_large$linewidth, na.rm = TRUE), mean(edge_small$linewidth, na.rm = TRUE))
    expect_gte(length(unique(edge_large$linetype)), 2)
    expect_equal(length(unique(edge_no_coupling$linetype)), 1)
})

test_that("cnetplot encodes mnsea feature sign with node colors", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "none")
    feature_nodes <- ggplot2::ggplot_build(p)$data[[2]]

    expect_gte(length(unique(feature_nodes$colour)), 2)
})

test_that("cnetplot distinguishes mnsea pathway and feature nodes with shapes", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "none")
    built <- ggplot2::ggplot_build(p)$data

    expect_true(all(built[[2]]$shape == 21))
    expect_true(all(built[[3]]$shape == 23))
})

test_that("cnetplot uses explicit mnsea legend titles", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "none")

    expect_equal(p$scales$get_scales("linetype")$name, "Edge type")
    expect_equal(p$scales$get_scales("shape")$name, "Node type")
    expect_equal(p$scales$get_scales("fill")$name, "Layer")
    expect_equal(p$scales$get_scales("colour")$name, "Feature sign")
    expect_equal(p$scales$get_scales("size")$name, "Feature magnitude")
})


test_that("gsInfo dispatches to mnseaResult with layer", {
    x <- mock_mnsea_result()

    d_collapsed <- enrichplot:::gsInfo(x, "T1")
    d_rna <- enrichplot:::gsInfo(x, "T1", layer = "rna")
    d_protein <- enrichplot:::gsInfo(x, "T1", layer = "protein")

    expect_true(all(c("ID", "Description", "runningScore", "position", "layer") %in% colnames(d_collapsed)))
    expect_equal(unique(as.character(d_collapsed$layer)), "collapsed")
    expect_equal(unique(as.character(d_rna$layer)), "rna")
    expect_equal(unique(as.character(d_protein$layer)), "protein")
    expect_equal(as.numeric(d_rna$geneList), as.numeric(x@layer_scores$rna))
})

test_that("gseaplot2 works for mnseaResult with collapsed scores", {
    x <- mock_mnsea_result()

    p <- gseaplot2(x, geneSetID = "T1")

    expect_s3_class(p, "gglist")
    expect_true(all(c("ID", "Description", "runningScore", "position", "layer") %in% colnames(p[[1]]$data)))
    expect_equal(unique(as.character(p[[1]]$data$layer)), "collapsed")
})

test_that("gseaplot2 mnseaResult supports single-layer ranked scores", {
    x <- mock_mnsea_result()

    p <- gseaplot2(x, geneSetID = "T1", layer = "protein")

    expect_s3_class(p, "gglist")
    expect_equal(unique(as.character(p[[1]]$data$layer)), "protein")
    expect_equal(as.numeric(p[[1]]$data$geneList), as.numeric(x@layer_scores$protein))
})

test_that("gseaplot2 mnseaResult resolves numeric geneSetID", {
    x <- mock_mnsea_result()

    p_idx <- gseaplot2(x, geneSetID = 1, layer = "rna")
    p_id <- gseaplot2(x, geneSetID = "T1", layer = "rna")

    expect_equal(as.numeric(p_idx[[1]]$data$runningScore), as.numeric(p_id[[1]]$data$runningScore))
    expect_equal(as.character(p_idx[[1]]$data$Description), as.character(p_id[[1]]$data$Description))
})

test_that("gsearank works for nsea and mnsea results", {
    x <- mock_mnsea_result()

    p_collapsed <- gsearank(x, geneSetID = "T1")
    p_rna <- gsearank(x, geneSetID = "T1", layer = "rna")
    tbl <- gsearank(x, geneSetID = "T1", output = "table", layer = "rna")

    expect_s3_class(p_collapsed, "ggplot")
    expect_s3_class(p_rna, "ggplot")
    expect_true(all(c("gene", "rank in geneList", "running ES", "core enrichment") %in% colnames(tbl)))
    expect_equal(unique(as.character(p_rna$data$layer)), "rna")
})

test_that("gseaplot2 multi-geneSetID works for mnsea collapsed scores", {
    x <- mock_mnsea_result()

    p <- gseaplot2(x, geneSetID = c("T1", "T2"))

    expect_s3_class(p, "gglist")
    expect_setequal(unique(as.character(p[[1]]$data$Description)), c("Pathway 1", "Pathway 2"))
})


test_that("pairwise_termsim works for mnseaResult", {
    x <- mock_mnsea_result()

    y <- pairwise_termsim(x, showCategory = 2)

    expect_s4_class(y, "mnseaResult")
    expect_true(is.matrix(y@termsim) || is.data.frame(y@termsim))
    expect_equal(nrow(y@termsim), 2)
    expect_equal(ncol(y@termsim), 2)
    expect_equal(y@method, "JC")
})

test_that("treeplot works for mnseaResult after pairwise_termsim", {
    x <- mock_mnsea_result()
    y <- pairwise_termsim(x, showCategory = 2)

    p <- treeplot(y, showCategory = 2)

    expect_s3_class(p, "ggtree")
})

test_that("treeplot supports single-pathway boundary cases", {
    x <- mock_mnsea_result()
    # Only one pathway with precomputed termsim should not error
    x@termsim <- matrix(1, nrow = 1, ncol = 1,
        dimnames = list("Pathway 1", "Pathway 1"))
    x@result <- x@result[1, , drop = FALSE]

    p <- treeplot(x, showCategory = 1)
    expect_s3_class(p, "ggplot")
})

test_that("barplot works for mnseaResult without falling back to graphics::barplot", {
    x <- mock_mnsea_result()

    p <- barplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
    expect_true(all(c("Description", "Count", "p.adjust") %in% colnames(p$data)))
})


test_that("hplot works for mnseaResult with collapsed scores", {
    x <- mock_mnsea_result()

    p <- hplot(x, geneSetID = "T1")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("ID", "Description", "runningScore", "position", "layer") %in% colnames(p$data)))
    expect_equal(unique(as.character(p$data$layer)), "collapsed")
})

test_that("hplot supports single-layer ranked scores", {
    x <- mock_mnsea_result()

    p <- hplot(x, geneSetID = "T1", layer = "protein")

    expect_s3_class(p, "ggplot")
    expect_equal(unique(as.character(p$data$layer)), "protein")
    expect_equal(as.numeric(p$data$geneList), as.numeric(x@layer_scores$protein))
})
