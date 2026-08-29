test_that("compute_rewiring_score returns NA for a single nseaResult without reference", {
    # mnsea mock represents the no-reference case as well
    x <- mock_mnsea_result()

    df <- enrichplot:::compute_rewiring_score(x)

    expect_s3_class(df, "data.frame")
    expect_true(all(c("ID", "Description", "layer", "rewiring_score", "leading_edge_overlap") %in% colnames(df)))
    expect_true(all(is.na(df$rewiring_score)))
})

test_that("compute_rewiring_score compares mnsea layers with a reference", {
    x <- mock_mnsea_result()

    df <- enrichplot:::compute_rewiring_score(
        x,
        selected_layer = "rna",
        reference_layer = "protein"
    )

    expect_true(all(df$leading_edge_overlap >= 0 & df$leading_edge_overlap <= 1))
    expect_equal(df$leading_edge_overlap, c(1, 1))
    expect_equal(df$rewiring_score, c(0, 0))
})

test_that("classify_mechanism_state returns deterministic labels", {
    state <- enrichplot:::classify_mechanism_state(
        nes_shift = c(0.1, 0.3, 0.5, 0.9),
        rewiring_score = c(0.2, 0.4, 0.6, 0.8)
    )

    expect_equal(state, c("conserved", "context_specific", "rewired", "rewired"))
    expect_true(all(is.na(enrichplot:::classify_mechanism_state(NA_real_, NA_real_))))
})

test_that("summarize_nsea_mechanism returns required columns", {
    x <- mock_mnsea_result()

    df <- enrichplot:::summarize_nsea_mechanism(x)

    expect_s3_class(df, "data.frame")
    expect_true(all(c(
        "ID", "Description", "NES", "p.adjust",
        "leading_edge_size", "leading_edge_overlap",
        "rewiring_score", "centrality_shift", "mechanism_class"
    ) %in% colnames(df)))
    expect_true(all(c("reference_NES", "delta_NES") %in% colnames(df)))
    expect_equal(df$leading_edge_size, c(2L, 2L))
    expect_true(all(is.na(df$rewiring_score)))
})

test_that("summarize_nsea_mechanism computes delta_NES from a reference result", {
    x <- mock_mnsea_result()
    # Create a modified result with different NES values.
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    df <- enrichplot:::summarize_nsea_mechanism(
        x2,
        reference = x,
        selected_layer = "rna",
        reference_layer = "protein"
    )

    expect_equal(df$delta_NES, c(0.4, 0.4))
    expect_equal(df$reference_NES, c(1.4, -1.1))
    # With rewiring = 0 and delta_NES = 0.4 (> nes threshold), the class is
    # context_specific rather than conserved, proving classify uses real NES shift.
    expect_equal(df$mechanism_class, c("context_specific", "context_specific"))
    expect_equal(df$rewiring_score, c(0, 0))
})

test_that("extract_rewiring_features requires reference_layer", {
    x <- mock_mnsea_result()

    expect_error(
        enrichplot:::extract_rewiring_features(
            x,
            pathway_id = "T1",
            selected_layer = "rna",
            reference_layer = NULL
        ),
        "reference_layer must be supplied"
    )
})

test_that("extract_rewiring_features returns shared/gained/lost status", {
    x <- mock_mnsea_result()

    df <- enrichplot:::extract_rewiring_features(
        x,
        pathway_id = "T1",
        selected_layer = "rna",
        reference_layer = "protein",
        score_change_threshold = 0
    )

    expect_true(all(c("Feature", "score", "abs_score", "sign", "status") %in% colnames(df)))
    expect_true(all(df$status %in% c("shared", "gained", "lost", "shifted")))
})

test_that("summarize_nsea_mechanism respects custom thresholds", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    # With default thresholds delta_NES = 0.4 is "context_specific".
    df_default <- enrichplot:::summarize_nsea_mechanism(
        x2,
        reference = x,
        selected_layer = "rna",
        reference_layer = "protein"
    )
    expect_equal(df_default$mechanism_class, c("context_specific", "context_specific"))

    # With a high nes threshold, the same delta_NES becomes "conserved".
    df_strict <- enrichplot:::summarize_nsea_mechanism(
        x2,
        reference = x,
        selected_layer = "rna",
        reference_layer = "protein",
        thresholds = list(rewire = 0.5, nes = 1)
    )
    expect_equal(df_strict$mechanism_class, c("conserved", "conserved"))
})

