test_that("phaseplot works with a real nseaResult", {
    x <- mock_nsea_result()

    p <- phaseplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
    expect_true(all(c("NES", "rewiring_score", "leading_edge_size", "mechanism_class") %in% colnames(p$data)))
})

test_that("consensusmap rejects a bare nseaResult", {
    x <- mock_nsea_result()

    expect_error(
        consensusmap(x),
        "Pass a named list of results"
    )
})

test_that("mechanismflow rejects a bare nseaResult", {
    x <- mock_nsea_result()

    expect_error(
        mechanismflow(x),
        "requires a named list of results"
    )
})

test_that("gseaplot2 works with nseaResult", {
    x <- mock_nsea_result()

    p <- gseaplot2(x, geneSetID = "T1")

    expect_s3_class(p, "gglist")
    expect_true(all(c("Description", "runningScore", "position") %in% colnames(p[[1]]$data)))
})

test_that("gsearank works with nseaResult", {
    x <- mock_nsea_result()

    p <- gsearank(x, geneSetID = "T1")

    expect_s3_class(p, "ggplot")
})

test_that("hplot works with nseaResult", {
    x <- mock_nsea_result()

    p <- hplot(x, geneSetID = "T1")

    expect_s3_class(p, "ggplot")
})

test_that("barplot works with nseaResult", {
    x <- mock_nsea_result()

    p <- barplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
    expect_true(all(c("Description", "Count", "p.adjust") %in% colnames(p$data)))
})
