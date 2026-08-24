test_that("phaseplot works for nseaResult", {
    # use mnsea mock to exercise the shared internal path
    x <- mock_mnsea_result()

    p <- phaseplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
    expect_true(all(c("NES", "rewiring_score", "leading_edge_size", "mechanism_class") %in% colnames(p$data)))
})

test_that("phaseplot works for mnseaResult", {
    x <- mock_mnsea_result()

    p <- phaseplot(x, showCategory = 2, selected_layer = "rna", reference_layer = "protein")

    expect_s3_class(p, "ggplot")
})

test_that("phaseplot handles empty result boundaries", {
    x <- mock_mnsea_result()
    x@result <- x@result[0, , drop = FALSE]

    expect_error(
        phaseplot(x, showCategory = 2),
        "No mnsea pathways available"
    )
})
