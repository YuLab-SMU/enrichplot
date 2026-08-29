loadNamespace("DOSE")
loadNamespace("clusterProfiler")
loadNamespace("enrichit")

mock_enrich_result <- function() {
    result <- data.frame(
        ID = c("T1", "T2"),
        Description = c("dup", "dup"),
        GeneRatio = c("1/2", "1/2"),
        BgRatio = c("1/10", "1/10"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.02, 0.03),
        qvalue = c(0.02, 0.03),
        geneID = c("g1/g2", "g2/g3"),
        Count = c(2L, 2L),
        stringsAsFactors = FALSE
    )
    rownames(result) <- result$ID

    methods::new(
        "enrichResult",
        result = result,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        organism = "mock",
        ontology = "mock",
        gene = c("g1", "g2", "g3"),
        keytype = "UNKNOWN",
        universe = character(),
        gene2Symbol = character(),
        geneSets = list(
            T1 = c("g1", "g2"),
            T2 = c("g2", "g3")
        ),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

mock_foldchange <- function() {
    c(g1 = 1, g2 = 2, g3 = 3)
}

mock_comparecluster_result <- function() {
    gc <- list(
        A = c("1", "2", "3"),
        B = c("2", "3", "4")
    )
    df <- data.frame(
        Cluster = c("A", "A", "B", "B"),
        ID = c("T1", "T2", "T1", "T2"),
        Description = c("dup", "other", "dup", "other"),
        geneID = c("1/2", "2/3", "2/3", "3/4"),
        Count = c(2L, 2L, 2L, 2L),
        p.adjust = c(0.01, 0.02, 0.03, 0.04),
        stringsAsFactors = FALSE
    )

    methods::new(
        "compareClusterResult",
        compareClusterResult = df,
        geneClusters = gc,
        fun = "mock",
        gene2Symbol = character(),
        keytype = "UNKNOWN",
        readable = FALSE,
        .call = quote(mock()),
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list(),
        organism = "mock"
    )
}

mock_mnsea_result <- function() {
    result <- data.frame(
        ID = c("T1", "T2"),
        Description = c("Pathway 1", "Pathway 2"),
        setSize = c(2L, 2L),
        enrichmentScore = c(0.8, -0.5),
        NES = c(1.4, -1.1),
        pvalue = c(0.01, 0.03),
        p.adjust = c(0.02, 0.04),
        qvalue = c(0.02, 0.04),
        rank = c(1L, 2L),
        leading_edge = c("tags=50%, list=40%, signal=30%", "tags=50%, list=20%, signal=10%"),
        core_enrichment = c("g1/g2", "g2/g3"),
        stringsAsFactors = FALSE
    )
    rownames(result) <- result$ID

    pathway_contribution <- data.frame(
        ID = c("T1", "T1", "T2", "T2"),
        Description = c("Pathway 1", "Pathway 1", "Pathway 2", "Pathway 2"),
        layer = c("rna", "protein", "rna", "protein"),
        contribution = c(0.7, 0.3, 0.4, 0.6),
        share = c(0.6, 0.4, 0.45, 0.55),
        n_feature = c(2L, 2L, 2L, 2L),
        stringsAsFactors = FALSE
    )

    feature_contribution <- data.frame(
        ID = c("T1", "T1", "T1", "T1", "T2", "T2", "T2", "T2"),
        Description = c(
            "Pathway 1", "Pathway 1", "Pathway 1", "Pathway 1",
            "Pathway 2", "Pathway 2", "Pathway 2", "Pathway 2"
        ),
        Feature = c("g1", "g2", "g1", "g2", "g2", "g3", "g2", "g3"),
        layer = c("rna", "rna", "protein", "protein", "rna", "rna", "protein", "protein"),
        score = c(1.2, 0.8, -0.4, 1.0, -0.7, -0.5, -0.3, -0.9),
        abs_score = c(1.2, 0.8, 0.4, 1.0, 0.7, 0.5, 0.3, 0.9),
        is_core = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
        stringsAsFactors = FALSE
    )

    coupling_table <- data.frame(
        from_layer = c("rna", "rna"),
        from_id = c("g1", "g2"),
        to_layer = c("protein", "protein"),
        to_id = c("g1", "g2"),
        weight = c(0.8, 0.6),
        stringsAsFactors = FALSE
    )

    rna_adj <- Matrix::Matrix(
        matrix(
            c(
                0, 1, 0,
                1, 0, 1,
                0, 1, 0
            ),
            nrow = 3,
            byrow = TRUE,
            dimnames = list(c("g1", "g2", "g3"), c("g1", "g2", "g3"))
        ),
        sparse = TRUE
    )
    protein_adj <- Matrix::Matrix(
        matrix(
            c(
                0, 1, 0,
                1, 0, 1,
                0, 1, 0
            ),
            nrow = 3,
            byrow = TRUE,
            dimnames = list(c("g1", "g2", "g3"), c("g1", "g2", "g3"))
        ),
        sparse = TRUE
    )

    methods::new(
        "mnseaResult",
        result = result,
        organism = "mock",
        setType = "mock",
        geneSets = list(
            T1 = c("g1", "g2"),
            T2 = c("g2", "g3")
        ),
        geneList = c(g1 = 1.2, g2 = 0.8, g3 = -0.5),
        keytype = "UNKNOWN",
        permScores = matrix(0, 0, 0),
        params = list(),
        gene2Symbol = character(),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list(),
        multilayer_network = list(
            intra_matrices = list(
                rna = rna_adj,
                protein = protein_adj
            )
        ),
        layer_scores = list(
            rna = c(g1 = 1.2, g2 = 0.8, g3 = -0.5),
            protein = c(g1 = -0.4, g2 = 1.0, g3 = -0.9)
        ),
        collapsed_scores = c(g1 = 0.8, g2 = 0.9, g3 = -0.7),
        layer_weights = c(rna = 1, protein = 1.5),
        coupling_table = coupling_table,
        mode = "signed",
        iterations = as.integer(12),
        restart_prob = 0.7,
        collapse_method = "weighted_mean",
        target_layer = "",
        output_space = "union",
        pathway_contribution = pathway_contribution,
        feature_contribution = feature_contribution
    )
}

mock_nsea_result <- function() {
    result <- data.frame(
        ID = c("T1", "T2"),
        Description = c("Network Pathway 1", "Network Pathway 2"),
        setSize = c(2L, 2L),
        enrichmentScore = c(0.8, -0.5),
        NES = c(1.4, -1.1),
        pvalue = c(0.01, 0.03),
        p.adjust = c(0.02, 0.04),
        qvalue = c(0.02, 0.04),
        rank = c(1L, 2L),
        leading_edge = c("tags=50%, list=40%, signal=30%", "tags=50%, list=20%, signal=10%"),
        core_enrichment = c("g1/g2", "g2/g3"),
        stringsAsFactors = FALSE
    )
    rownames(result) <- result$ID

    adj <- Matrix::Matrix(
        matrix(
            c(
                0, 1, 0,
                1, 0, 1,
                0, 1, 0
            ),
            nrow = 3,
            byrow = TRUE,
            dimnames = list(c("g1", "g2", "g3"), c("g1", "g2", "g3"))
        ),
        sparse = TRUE
    )

    methods::new(
        "nseaResult",
        result = result,
        organism = "mock",
        setType = "mock",
        geneSets = list(
            T1 = c("g1", "g2"),
            T2 = c("g2", "g3")
        ),
        geneList = c(g1 = 1.2, g2 = 0.8, g3 = -0.5),
        keytype = "UNKNOWN",
        permScores = matrix(0, 0, 0),
        params = list(),
        gene2Symbol = character(),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list(),
        network = adj,
        diffusion_scores = c(g1 = 1.2, g2 = 0.8, g3 = -0.5),
        mode = "signed",
        iterations = as.integer(15),
        restart_prob = 0.7
    )
}
