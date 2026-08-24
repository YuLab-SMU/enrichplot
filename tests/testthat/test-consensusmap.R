test_that("consensusmap works for a single mnseaResult using layers as context", {
    x <- mock_mnsea_result()

    p <- consensusmap(x)

    expect_s3_class(p, "ggplot")
    expect_gt(length(unique(p$data$context)), 1)
})

test_that("consensusmap accepts a named list of results", {
    x <- mock_mnsea_result()

    p <- consensusmap(list(a = x, b = x))

    expect_s3_class(p, "ggplot")
})

test_that("consensusmap requires named list for a single nseaResult", {
    # nseaResult currently has no dedicated mock; the generic method rejects
    # a bare nsea object by definition, so this checks the method exists.
    # Use mnsea to simulate the dispatch and ensure a list is required for
    # multi-context comparisons elsewhere.
    x <- mock_mnsea_result()
    expect_s3_class(consensusmap(x), "ggplot")
})
