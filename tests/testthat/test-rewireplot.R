test_that("rewireplot works for mnseaResult", {
    x <- mock_mnsea_result()

    p <- rewireplot(
        x,
        pathway_id = "T1",
        selected_layer = "rna",
        reference_layer = "protein"
    )

    expect_s3_class(p, "ggplot")
})

test_that("rewireplot requires reference_layer", {
    x <- mock_mnsea_result()

    expect_error(
        rewireplot(
            x,
            pathway_id = "T1",
            selected_layer = "rna",
            reference_layer = NULL
        ),
        "reference_layer must be supplied"
    )
})

test_that("rewireplot requires selected_layer", {
    x <- mock_mnsea_result()

    expect_error(
        rewireplot(
            x,
            pathway_id = "T1",
            selected_layer = NULL,
            reference_layer = "protein"
        ),
        "selected_layer must be supplied"
    )
})

test_that("rewireplot maps feature status from extract_rewiring_features", {
    x <- mock_mnsea_result()
    status <- enrichplot:::extract_rewiring_features(
        x,
        pathway_id = "T1",
        selected_layer = "rna",
        reference_layer = "protein",
        score_change_threshold = 0
    )

    expect_true(all(status$status %in% c("shared", "gained", "lost", "shifted")))
})
