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

test_that("phaseplot uses delta_NES when reference is supplied", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- phaseplot(
        x2,
        reference = x,
        selected_layer = "rna",
        reference_layer = "protein",
        showCategory = 2
    )

    expect_s3_class(p, "ggplot")
    expect_true(all(c("delta_NES", "rewiring_score", "leading_edge_size", "mechanism_class") %in% colnames(p$data)))
    expect_false(all(is.na(p$data$delta_NES)))
})

test_that("phaseplot can use legacy NES axis", {
    x <- mock_mnsea_result()

    p <- phaseplot(x, showCategory = 2, x_axis = "NES")

    expect_s3_class(p, "ggplot")
    expect_true("NES" %in% colnames(p$data))
})

test_that("phaseplot rejects delta_NES without reference", {
    x <- mock_mnsea_result()

    expect_error(
        phaseplot(x, showCategory = 2, x_axis = "delta_NES"),
        "No `delta_NES` available"
    )
})

test_that("phaseplot respects custom thresholds", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p_strict <- phaseplot(
        x2,
        reference = x,
        selected_layer = "rna",
        reference_layer = "protein",
        showCategory = 2,
        thresholds = list(rewire = 0.5, nes = 1)
    )

    expect_s3_class(p_strict, "ggplot")
    expect_true("mechanism_class" %in% colnames(p_strict$data))
    expect_equal(unique(p_strict$data$mechanism_class), "conserved")
})
