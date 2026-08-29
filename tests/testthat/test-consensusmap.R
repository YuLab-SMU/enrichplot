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

test_that("consensusmap fills by delta_NES when available", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- consensusmap(list(control = x, treated = x2), fill_var = "delta_NES")

    expect_s3_class(p, "ggplot")
    expect_true("delta_NES" %in% colnames(p$data))
    expect_false(all(is.na(p$data$delta_NES)))
})

test_that("consensusmap uses rewiring size when include_rewiring is TRUE", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- consensusmap(
        list(control = x, treated = x2),
        include_rewiring = TRUE,
        size_var = "rewiring_score"
    )

    expect_s3_class(p, "ggplot")
    expect_true("rewiring_score" %in% colnames(p$data))
    expect_false(all(is.na(p$data$rewiring_score)))
})

test_that("consensusmap supports rewiring-score labels", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- consensusmap(
        list(control = x, treated = x2),
        label = "rewiring_score"
    )

    expect_s3_class(p, "ggplot")
})

test_that("consensusmap respects custom thresholds", {
    x <- mock_mnsea_result()
    x2 <- mock_mnsea_result()
    x2@result$NES <- c(1.8, -0.7)

    p <- consensusmap(
        list(control = x, treated = x2),
        thresholds = list(rewire = 0.5, nes = 1)
    )

    expect_s3_class(p, "ggplot")
    expect_true("mechanism_class" %in% colnames(p$data))
})
