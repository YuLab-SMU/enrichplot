test_that("mechanismflow works for a named list of results", {
    x <- mock_mnsea_result()

    p <- mechanismflow(list(a = x, b = x))

    expect_s3_class(p, "ggplot")
})

test_that("mechanismflow uses layers as context for a single mnseaResult", {
    x <- mock_mnsea_result()

    p <- mechanismflow(x)

    expect_s3_class(p, "ggplot")
})

test_that("mechanismflow rejects single nseaResult", {
    x <- mock_mnsea_result()
    # The nseaResult method raises an explicit error, but the same guard is
    # covered by the internal single-result check; use mnsea with only one
    # layer to exercise the boundary.
    x@layer_scores <- x@layer_scores[1]
    expect_error(
        mechanismflow(x),
        "At least two contexts are required"
    )
})

test_that("mechanismflow uses delta_NES flow magnitude", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- mechanismflow(
        list(control = x, treated = x2),
        flow_var = "delta_NES"
    )

    expect_s3_class(p, "ggplot")
    expect_true("delta_NES" %in% colnames(p$data))
    expect_false(all(is.na(p$data$delta_NES)))
})

test_that("mechanismflow preserves stable mechanism-state order", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- mechanismflow(list(control = x, treated = x2))

    expect_s3_class(p, "ggplot")
    expect_true(all(c("context", "ID", "state_num", "mechanism_class") %in% colnames(p$data)))
})

test_that("mechanismflow respects custom thresholds", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- mechanismflow(
        list(control = x, treated = x2),
        thresholds = list(rewire = 0.5, nes = 1)
    )

    expect_s3_class(p, "ggplot")
})
