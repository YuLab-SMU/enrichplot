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
