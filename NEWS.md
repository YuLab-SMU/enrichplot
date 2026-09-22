# enrichplot 1.99.9

+ fix `treeplot()` split-aware faceting for GSEA results (#169): `split` is now carried into tree, tip, and clade metadata so `treeplot(..., split = ".sign") + facet_grid(. ~ .sign)` builds instead of dropping the faceting variable from every layer (2026-09-22, Tue)

# enrichplot 1.99.8

+ fix `treeplot()` cluster color assignment when `nCluster` reaches two digits (#171): cluster palettes and clade-label groups now follow numeric cluster ids instead of lexical ordering, so groups like `cluster_10` no longer steal `cluster_2` colors when the cluster count increases (2026-09-22, Tue)

# enrichplot 1.99.7

+ fix `emapplot()` compare-cluster pie nodes across ontologies (#228): ontology-specific terms that share the same `Description` now keep stable ID-backed labels all the way into pie-layer data alignment, so `compareCluster(..., ont = "ALL")` no longer collapses those nodes or breaks while building the pie overlay (2026-09-22, Tue)

# enrichplot 1.99.6

+ fix `goplot()` DAG construction for top-level GO terms: parent edges that point to the synthetic `all` root are now dropped before the graph is built, so plots that include terms such as `GO:0008150` no longer fail with `Some vertex names in \`d\` are not listed in \`vertices\`` (2026-09-22, Tue)

# enrichplot 1.99.5

+ fix `heatplot()` dot-mode p-value scaling: zero or non-positive gene p-values are now clamped to the smallest positive double before the reversed log-size transform is applied, so significance-sized dots no longer emit infinite-value warnings for exact-zero inputs (2026-09-22, Tue)
+ fix `ridgeplot()` blank rows for undersized core gene sets (#288): pathways with fewer than three ranked values are now dropped before `geom_density_ridges()` is built, and the function now errors clearly when no selected pathway has enough values to estimate a density, so two-gene core sets no longer leave empty y-axis slots in the plot (2026-09-21, Mon)

# enrichplot 1.99.4

+ fix `cnetplot()` for `compareClusterResult` terms with duplicated descriptions (#279): category nodes now use stable ID-backed labels internally, so distinct terms that share the same `Description` are no longer merged into one network node, with regression coverage for the duplicated-label case (2026-09-21, Mon)

# enrichplot 1.99.3

+ fix `treeplot()` heatmap panels for `compareClusterResult`: the `cluster_panel = "heatMap"` path now calls `ggtree::gheatmap()` with the active tree plot object instead of treating it like a regular layer, so compare-cluster treeplots render again instead of failing with a missing `data` argument error (2026-09-21, Mon)
+ fix `treeplot()` dotplot panels for `compareClusterResult` (#232, #224): the `cluster_panel = "dotplot"` path now passes `ggtreeExtra::geom_fruit()` the plain term columns it expects, so compare-cluster treeplots render again on current `ggtreeExtra` builds instead of failing while decoding the `y` mapping (2026-09-21, Mon)

# enrichplot 1.99.2

+ support importing results from external enrichment tools: `import_enrichr()`, `import_gprofiler2()`, `import_webgestalt()` and `import_fgsea()` map enrichr / g:Profiler / WebGestaltR / fgsea output tables to `enrichResult` / `gseaResult` objects that plug into the 'enrichplot' visualization functions; the 'enrichit' constructors `as_enrichResult()` / `as_gseaResult()` are re-exported for other table formats (2026-09-21, Mon)
+ fix `dotplot()` legend keys under plot composition (#273): size legends now keep the hollow point shape after `cowplot::plot_grid()` / similar grob composition, instead of reverting to solid circles in combined figures (2026-09-21, Mon)
+ fix `dotplot()` size scaling for enrichment results (#118): `size = "Percentage"` now derives a percentage column from `GeneRatio` for `enrichResult` / `gseaResult` data instead of failing at draw time with a missing-column error (2026-09-21, Mon)
+ fix `barplot()` width handling (#201): `width = ...` is now forwarded to the internal `geom_col()` layer for both enrichment-result and compare-cluster barplots, so bar thickness can be adjusted directly without stacking a second `geom_col()` layer on top of the original bars (2026-09-21, Mon)
+ fix `dotplot()` selection ordering for numeric `showCategory` (#345, #219): the function now orders the fortified data by `orderBy` first and only then takes the requested top rows, so the leading categories stay stable when `showCategory` changes and `orderBy` is honored correctly (2026-09-21, Mon)

# enrichplot 1.99.1

+ fix `upsetplot()` for readable `gseaResult` objects (#179): the fold-change vector is now remapped through `fc_readable()`, so `setReadable()` results no longer lose all ranked-score values when pathway genes are shown as symbols (2026-09-21, Mon)
+ improve `gseaplot()` / `gseaplot2()` multi-panel compatibility with `cowplot` (#239): the `gglist`-level `cowplot::as_grob()` bridge is now implemented in `aplot`, so `plot_grid()` / `ggarrange()` work when paired with an `aplot` version that provides that helper, without making it a hard requirement for `enrichplot` itself (2026-09-21, Mon)
+ fix `upsetplot()` boxplot overlays for `gseaResult` and `mnseaResult` (#178): the boxplot layer now suppresses its own outlier glyphs so jittered feature points are drawn only once instead of being duplicated on top of boxplot outliers (2026-09-21, Mon)
+ fix `pairwise_termsim()` for `enrichResult` objects whose raw result table has terms but the object cutoffs filter them all out of `as.data.frame()` (#269): term selection now uses the raw result rows, so `showCategory` can still pick the requested top terms and downstream plots such as `emapplot()` continue to work for non-significant result tables (2026-09-21, Mon)
+ fix grouped `emapplot()` / `ssplot()` legend control (#292): the compatibility arguments `group` and `group_legend` are accepted again, and grouped layouts no longer force the "groups" legend on when `group_legend = FALSE` is requested (2026-09-21, Mon)
+ fix `gseaplot2()` hit-bin rectangles for single gene sets (#221, #20): the colored bins under the hit ticks now follow the ranked-list direction instead of mirroring the cumulative hit counts, so highly one-sided enrichments no longer collapse the wide interval onto the wrong end of the plot (2026-09-21, Mon)
+ fix `cnetplot()` / `emapplot()` for `compareClusterResult` pie nodes when duplicated `(Cluster, Description)` rows are present (#314): pie counts are now aggregated before widening, avoiding list-columns and the tidyr cast error ("Can't convert `fill` <double> to <list>"), with regression coverage for duplicated cluster-term inputs (2026-09-21, Mon)
+ fix `barplot()` for `compareClusterResult` objects: the default `by = "geneRatio"` and `by = "rowPercentage"` crashed in `plotting.clusterProfile()`, and `by = "count"` failed at rendering time under `ggplot2` 4.x; `by` is now mapped to the fortify-produced column and bars are drawn with `geom_col()` (2026-09-21, Mon)
+ fix `emapplot()` / `ssplot()` with similarity measures other than 'JC' (e.g., 'Wang') (#309): label-keyed similarity matrices were re-mapped as term IDs, producing NA edges ("edge data frame contains NAs"); the stale re-mapping in `build_emap_graph()` was removed (2026-09-21, Mon)
+ add a plotting regression suite (`test-plotting-regression.R`) covering the tutorial-facing visualization functions, including a dispatch canary for the `ggplot() + theme_dose()` failure seen under `ggplot2` 4.0.x with S7 < 0.2.2; ggplot outputs are evaluated with `ggplot_build()` to catch bad aesthetics (2026-09-21, Mon)

# enrichplot 1.33.1

+ complete remaining mechanism-plot enhancements: `pairwise_termsim()` now supports layer-aware similarity for `mnseaResult`, classification thresholds are exposed through `phaseplot()` / `consensusmap()` / `mechanismflow()`, and an explicit `nseaResult` mock plus coverage has been added for nsea plotting paths; `gsInfo.gseaResult()` also defaults `exponent` to 1 when `params` lacks it (2026-08-29, Sun)
+ refactor `gsInfo()` into an S3 generic and add `layer`-aware running-score support to `gseaplot2()`, `gsearank()` and `hplot()` for `nseaResult` / `mnseaResult`; `hplot()` is now implemented with base `ggplot2` geoms and no longer requires `ggHoriPlot` (2026-08-23, Sun)
+ add `pairwise_termsim()` support for `mnseaResult` so `treeplot()`, `emapplot()` and `ssplot()` share one layer-aware similarity definition, with single-pathway treeplot boundary handling
+ add `barplot.gseaResult()` so `nseaResult` / `mnseaResult` no longer fall through to `graphics::barplot.default`
+ add mechanism-oriented helper layer (`compute_rewiring_score()`, `classify_mechanism_state()`, `summarize_nsea_mechanism()`, `extract_rewiring_features()`) with deterministic tests
+ add `phaseplot()` for enrichment-shift versus rewiring overviews and `rewireplot()` for pathway-specific feature-level rewiring evidence
+ add `consensusmap()` for multi-context mechanism agreement and `mechanismflow()` for pathway state transitions across layers/conditions
+ refine mechanism plots with real cross-object comparisons: `summarize_nsea_mechanism()` now accepts a `reference` result to compute `delta_NES` and cross-object `rewiring_score`; `phaseplot()` supports `reference` / `x_axis` / `size_var`; `consensusmap()` uses fill for NES/delta NES and point size for rewiring/overlap; `mechanismflow()` uses flow magnitude and stable mechanism-state ordering (2026-08-23, Sun)
+ add a minimal `ssplot.mnseaResult()` that projects selected pathways into a similarity-space overview using layer-aware feature overlap, while reusing `emapplot()` semantics and adding stable fallbacks for one- or two-pathway layouts (2026-06-25, Thu)
+ add a minimal `upsetplot.mnseaResult()` that summarizes shared feature overlaps across selected pathways with collapsed-score or single-layer views, including support for score magnitude display and `core_enrichment` filtering (2026-06-25, Thu)
+ add a minimal `ridgeplot.mnseaResult()` that shows pathway-level feature score distributions from collapsed scores or a selected single layer, with regression coverage for layer-aware ranked scores and `core_enrichment` filtering (2026-06-25, Thu)
+ add a minimal `gseaplot.mnseaResult()` that supports collapsed-score and single-layer running-score views for one pathway at a time, with regression coverage for stable pathway selection and layer-aware ranked scores (2026-06-25, Thu)
+ batch-refine `mnsea` plot semantics by aligning `layer` filtering and readable legend labels across `dotplot()`, `heatplot()`, `cnetplot()` and `emapplot()`, while fixing `emapplot.mnseaResult()` to retain all selected pathways when rebuilding overlap graphs after layer filtering, with expanded regression coverage (2026-06-25, Thu)
+ add a minimal `emapplot.mnseaResult()` that reuses cached term similarity when available and otherwise falls back to internal `JC` overlap for pathway-level map plots, with regression coverage (2026-06-24, Wed)
+ batch-refine `cnetplot.mnseaResult()` readability by splitting pathway and feature label layers, preferring shared features when labels are capped, and stabilizing layer ordering with expanded regression coverage (2026-06-24, Wed)
+ refine default label selection in `cnetplot.mnseaResult()` to keep pathway annotations while deduplicating repeated feature labels across layers, with regression coverage for the quieter defaults (2026-06-24, Wed)
+ clarify `cnetplot.mnseaResult()` legend titles for edge type, node type, layer, feature sign, and feature magnitude, with regression coverage for the updated defaults (2026-06-24, Wed)
+ distinguish pathway and feature nodes in `cnetplot.mnseaResult()` with explicit node-type shapes and regression coverage for the updated legend semantics (2026-06-24, Wed)
+ refine `cnetplot.mnseaResult()` with edge-type legends, effective `size_edge` scaling, and feature-node sign encoding backed by lightweight regression tests (2026-06-24, Wed)
+ align default `pathway_id` resolution across `mnsea` helpers and feature-level `heatplot()`, and add `share` / `exclusive` label support to `cnetplot.mnseaResult()` with regression coverage (2026-06-24, Wed)
+ add `cnetplot.mnseaResult()` for pathway-specific multilayer subnetworks, including pathway anchor nodes and lightweight regression coverage for the new network view (2026-06-24, Wed)
+ add `heatplot.mnseaResult()` for term-layer and pathway-specific feature heatmaps, and cover the new `mnsea` helper/plotting workflow with lightweight tests (2026-06-24, Wed)
+ refactor shared plot data preparation for `cnetplot()`, `emapplot()`, `heatplot()` and `pairwise_termsim()` around unified term selection helpers, and add smoke tests for `compareClusterResult` network visualizations (2026-06-24, Wed)
+ add a minimal `testthat` skeleton for regression coverage, and align `update_n()` / `pairwise_termsim()` / `get_similarity_matrix()` with stable term selection semantics (2026-06-24, Wed)
+ fix `heatplot(showTop)` to fail early when `foldChange` is missing, correct the `reverse` behavior in `set_enrichplot_color()`, and add runtime checks for optional plotting dependencies (2026-06-24, Wed)
+ harden term selection and label handling across `cnetplot()`, `emapplot()`, `pairwise_termsim()` and `upsetplot()` by using stable term identifiers internally while keeping display labels readable (2026-06-24, Wed)

# enrichplot 1.32.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# enrichplot 1.31.5

+ `cnetplot.compareClusterResult()` now supports `categorySizeBy` for category pie sizing and aligns docs with `ggtangle::cnetplot()` semantics (2026-04-22, Wed)
+ `ridgeplot` now supports `stat` parameter (default is 'density_ridges' and can be changed to 'binline') (2026-04-01, Wed, #343)
+ manhattan plot for enriched result (2026-03-26, Thu)
+ update roxygen document to use markdown syntax (2026-03-02, Mon)
+ bug fixed in xy-lab format in `ssplot()` (2026-03-02, Mon)
+ bug fixed in formula supports in `dotplot()` (2026-02-26, Thu)

# enrichplot 1.31.4

+ fix `cnetplot()` S3 generic/method consistency warnings (2026-01-14, Wed)
+ fix `treeplot()` column selection bug when color variable equals size variable (2026-01-14, Wed)
+ fix `fortify.compareClusterResult()` warnings about missing imports and global variables (2026-01-14, Wed)
+ remove `plyr` and use `dplyr` in `method-fortify.R` (2026-01-14, Wed)
+ fixed `treeplot()` issue where `pairwise_termsim()` with method="JC" produced unnamed similarity matrix, causing "undefined column selected" error (2025-01-14)
+ fixed `fortify.compareClusterResult()` warning "NAs introduced by coercion" when Cluster names are not numeric (2025-01-14)
+ bug fixed in `barplot()` as `fortify()` generic in `ggplot2` checks for unused arguments in `...` (2026-01-14, Wed)
+ remove `categorySize` parameter in `cnetplot()` (2026-01-14, Wed)
+ bug fixed in `goplot()` as `GOSemSim` uses cache (2026-01-13, Tue)
  - also fix `gotbl` object not found issue (2026-01-13, Tue)
+ re-export `geneID`, `geneInCategory` and `gseaScores` from 'enrichit' (2026-01-12, Mon)
+ update documentation: fix typos, grammar errors and use modern markdown syntax (2026-01-12, Mon)
+ bug fixed in `update_n()` if `showCategory` is a vector of term names (2026-01-08, Thu)
+ avoid the "condition has length > 1" error in `outer()` by using `Vectorize()` (2026-01-08, Thu)

# enrichplot 1.31.3

+ use 'enrichit' package (2025-12-07, Sun)
+ optimize source code (2025-12-02, Tue)
+ error handling functions imported from 'yulab.utils' (2025-12-01, Mon)

# enrichplot 1.31.2

+ add 'fc_threshold' parameter to `cnetplot` (2025-11-30, Sun, #338)
  - requires 'ggtangle' v>= 0.0.9
+ update all line width aes mapping from 'size' to 'linewidth' (2025-11-30, Sun)
+ add 'node_label_size' parameter for `emapplot` (2025-11-30, Sun)
+ remove `emapplot` parameters, 'group', 'group_style' and 'label_group_style' (#339) 
+ add 'showTop' parameter to limit number of genes shown in `heatplot()` and distinguish tip point size variable for `treeplot()` through internal parameter `size_var` (2025-11-23, Sat, #335) 

# enrichplot 1.31.1

+ import `ggfun::%<+%` (2025-11-18, Tue)
+ update `ssplot()`, `treeplot()` and `get_wordcloud()` (2025-11-15, Sat)
+ change `set_enrichplot_color(transform = 'identity')` as default behavior (2025-11-11, Tue)
  - now it only sets the color scale without changing the transform method
  - explicitly set `transform = 'log10'` in `dotplot`
+ use 'quarto' as vignette engine (2025-11-11, Tue)
+ use `set_enrichplot_color(transform = 'identity')` in `heatplot` (2025-11-11, Tue)
+ use `set_enrichplot_color(transform = 'identity')` in `cnetplot` (2025-11-05, Wed)

# enrichplot 1.30.0

+ Bioconductor RELEASE_3_22 (2025-11-01, Sat)

# enrichplot 1.29.4

+ remove deprecated `aes_string`/`aes_` (2025-10-23, Thu, #332)

# enrichplot 1.29.3

+ bug fixed of `cnetplot` for `CompareClusterResult` (2025-09-13, Sat, #329)
  - color gene according to the gene cluster info
+ bug fixed in pie scale label (2025-07-14, Mon, #328)

# enrichplot 1.29.2

+ update `treeplot` with two parameters, `leave_fontsize` and `clade_fontsize` (2025-07-12, Sat, #324, #325)
  - remove the `fontsize` parameter as it only works for `clade_fontsize`
+ 'log10' transform for pvalue color scale by default (2025-07-12, Sat, #316)
+ introduce new parameters in `gseaplot2()` (2025-07-12, Sat)
  - `pvalue_table_columns`
  - `pvalue_table_rownames`
  - <https://github.com/YuLab-SMU/clusterProfiler/issues/774>
  
# enrichplot 1.29.1

+ throw error in `goplot()` if ontology is not one of the 'MF', 'CC' or 'BP' (2025-04-28, Mon, clusterProfiler#768)

# enrichplot 1.28.0

+ Bioconductor RELEASE_3_21 (2025-04-17, Thu)

# enrichplot 1.27.5

+ able to scale pie size for 'compareClusterResult' (2025-03-11, Tue, #308, #311)

# enrichplot 1.27.4

+ adjust pie size and category label position in `cnetplot()` (2025-01-08, Wed, #306)
+ clean up code (2024-12-20, Fri)

# enrichplot 1.27.3

+ scale pies and add pie legend in `emapplot()` (2024-12-12, Thu, #304)
+ a safe way to extract gene sets in `ridgeplot()` (2024-12-12, Thu, #303)

# enrichplot 1.27.2

+ `emapplot()` now allows passing color to a specific color, e.g., color = "black" (2024-11-29, Fri, #300)
+ bug fixed in `emapplot()` 
  - `size_category` now works for pie node (2024-11-29, Fri, #301)
  - legend of term nodes will be retained when `group = TRUE` (2024-11-29, Fri, #300)
+ supports passing ID to 'showCategory' in `ridgeplot()` (2024-11-06, Wed, #295)
+ enhancement of `cnetplot()` (2024-11-06, Wed)
  - 'node_label' can be a vector of selected items/genes to specify the items to be displayed (#293)
  - 'node_label' can be 'exclusive' to label genes that are uniquely belongs to categories (#253)
  - 'node_label' can be 'share' to label genes that are share between categories (#253)
  - 'node_label' can be, e.g. '> 1' or '< 1', to label genes that have log2FC values larger or smaller than the threshold (#253) 
  - supports using `ggtangle::geom_cnet_label()` to label items/genes in independent layer (#194, #266, #267)
+ fixed `ridgeplot()` when selecting a specific gene set and plotting non-core genes (2024-11-06, Wed, #298)

# enrichplot 1.27.1

+ add 'ID' parameter in `goplot()` (2024-10-30, Wed)
  - <https://github.com/YuLab-SMU/enrichplot/issues/292#issuecomment-2445788948>

# enrichplot 1.26.0

+ Bioconductor RELEASE_3_20 (2024-10-30, Wed)

# enrichplot 1.25.6

+ pretty gene count legend (2024-10-29, Tue, #271)

# enrichplot 1.25.5

+ new `emaplot()`, `goplot()`, `cnetplot()` and `ssplot()`, all power by 'ggtangle' package (2024-10-24, Thu)
+ re-export `ggtangle::cnetplot()` (2024-10-24, Thu)
+ remove `drag_network()` (2024-10-24, Thu)

# enrichplot 1.25.4

+ fixed `goplot()` (2024-10-23, Wed, #297, #732, #718)

# enrichplot 1.25.3

+ `hplot()`: Horizontal plot for GSEA result (2024-08-27, Tue)

# enrichplot 1.25.2

+ fixed bug in `ridgeplot()` (2024-08-19, Mon, clusterProfiler#704)

# enrichplot 1.25.1

+ fixed GeneRatio in dotplot as character of fraction issue (2024-08-16, Fri, clusterProfiler#715)
+ use `yulab.utils::yulab_msg()` for startup message (2024-07-26, Fri)
+ `dotplot2` to compare two selected clusters in 'compareClusterResult' object (2024-06-15, Sat)
+ `volplot` to visualize ORA result using volcano plot (2024-06-13, Thu)

# enrichplot 1.24.0

+ Bioconductor RELEASE_3_19 (2024-05-15, Wed)

# enrichplot 1.23.2

+ separate the JC similarity method (2023-12-11, Mon, #265)
+ fix the issue in `ridgeplot(showCategory)` : support a vector of Description, not ID(2023-12-1, Fri, #193)

# enrichplot 1.23.1

+ `ridgeplot()` supports passing a vector of selected pathways via the 'showCategory' parameter (2023-11-30, Thu, #193)
+ fix `treeplot()` to compatible with the current version of ggtree and ggtreeExtra. (2023-10-28, Sat)
+ add clusterPanel.params[["colnames_angle"]] parameter to set the angle of colnames. (2023-10-28, Sat)

# enrichplot 1.22.0

+ Bioconductor RELEASE_3_18 (2023-10-25, Wed)

# enrichplot 1.21.3

+ `set_enrichplot_color()`, a helper function to set colors (2023-09-13, Wed)
  - change default color: from c("red", "blue") to c("#e06663", "#327eba")
+ use `check_installed()` to check package dependency (2023-09-08, Fri, #254)

# enrichplot 1.21.2

+ introduce 'facet' parameter in `dotplot()` method for `compareClusterResult`. If `facet = "intersect"`, the dots will be separated by enriched pathway intersection among clusters. It can set to other variable that can be used for splitting the figure (e.g., "category" for KEGG results) (2023-08-21, Mon)

# enrichplot 1.21.1

+ fixed `cnetplot.compareClusterResult()` for only contains one cluster (2023-05-24, Wed, #243)

# enrichplot 1.20.0

+ Bioconductor RELEASE_3_17 (2023-05-03, Wed)

# enrichplot 1.19.2

+ fix `emapplot()` for parameter mismatch (2023-02-20, Mon)
+ fix `ridgeplot` for error when x@readable == TRUE and length(x@gene2Symbol) = 0 (2022-12-5, Mon)
+ fix `ridgeplot` for error when `x@readable == TRUE` and `length(x@gene2Symbol) = 0` (2022-12-5, Mon, #217)

# enrichplot 1.19.1

+ fix `cnetplot()` for `node_label` parameter is flipped(2022-12-04, Sun, #216)
+ bug fixed in `treeplot()`  (2022-11-18, Fri) 
+ enable `dotplot()` and `autofacet()` for `gseaResultList` object

# enrichplot 1.18.0

+ Bioconductor RELEASE_3_16 (2022-11-02, Wed)

# enrichplot 1.17.4

+ rename parameters of `emapplot()`, `centplot()` and  `treeplot()` (2022-09-11, Sun)

# enrichplot 1.17.3

+ align the dots in `treeplot()` (2022-10-1, Sat)
+ fix a bug in color legend of `treeplot()` (2022-10-1, Sat)

# enrichplot 1.17.2

+ `autofacet` to automatically split `barplot` and `dotplot` into several facets (2022-09-06, Tue)
+ `dotplot` method for `enrichResultList` object 
+ add parameters `hilight_category`, `alpha_hilight`, `alpha_nohilight` for `cnetplot()` and `emapplot` (2022-09-4, Sun)
+ change round digits of cnetplot scatterpie legend to 1 (2022_8_29, Mon).
+ `gsearank()` can export result as a table when `output = "table"` (2022-08-29, Mon, #184)
+ fix a bug in `fc_readable()` (2022-08-29, Mon, #189)
+ allows passing `color="NES"` to `dotplot()` for `gseaResult` object (2022-08-29, Mon, #14)

# enrichplot 1.17.1

+ fix a bug in https://github.com/YuLab-SMU/clusterProfiler/issues/488 (2022-08-25, Thu)
+ support multiple gene sets in `geom_gsea_gene` layer (2022-08-25, Thu)
+ `geom_gsea_gene` layer (2022-08-24, Wed)
+ add parameters `symbol` and `pvalue` for `heatplot.enrichResult()` (2022-08-20, Sat)
+ change default values of `group_category` and `node_label` in `ssplot()` (2022-07-04, Mon)
+ update document of `ssplot()` (2022-07-04, Mon)
+ `gseaplot()` and `gseaplot2()` return `gglist` object instead of plotting the figure (2022-05-05, Thu)
+ fix `ridgeplot` when `x@readable = TRUE` (2022-04-30, Sat)

# enrichplot 1.16.0

+ Bioconductor 3.15 release

# enrichplot 1.15.4

+ update `treeplot`: support passing rel object to `offset` and `offset_tiplab` (2022-04-24, Sun)

# enrichplot 1.15.3

+ export `drag_network' (2022-03-07, Mon)
+ update `cnetplot.enrichResult` to be supported by `drag_network`(2022-3-6, Sun)
+ add function `drag_network` to drag the nodes of networks (2022-2-25, Fri)
+ fix a bug in `goplot`: `goplot.gseaResult` need `setType` slot instead of `ontology` slot (2022-2-22, Tue)
+ return `gg` object instead of print it in `dotplot.compareClusterResult()` (2022-01-05, Wed, @altairwei, #160)

# enrichplot 1.15.2

+ add `label_format_tiplab` and `label_format_cladelab` parameters for `treeplot`(2021-12-24, Fri)
+ support treeplot of compareCluster(GSEA algorithm) result(2021-12-13, Mon)
+ support visualization of compareCluster(GSEA algorithm) result(2021-12-11, Sat)
+ support scientific notation for `gseaplot2`(2021-12-4, Sat)

# enrichplot 1.15.1

+ fixed R check by importing `utils`

# enrichplot 1.14.0

+ Bioconductor 3.14 release

# enrichplot 1.13.2

+ mv `ep_str_wrap` to `yulab.utils::str_wrap` (2021-10-13, Wed) 
+ adjust the order of legends for `dotplot`, `emapplot`, `cnetplot` and `treeplot`(2021-10-8, Fri)
+ update `treeplot`: add "dotplot" and "heatmap" panels for `treeplot`(2021-9-15, Wed)
+ update `dotplot`: enable `size` parameter applicable to other columns of compareClusterResult(2021-9-17, Fri)
+ enable `label_format` parameter for `heatplot` (2021-09-01, Wed)
+ add `get_ggrepel_segsize` function to set `segment.size` value for `ggrepel`(2021-08-29, Sun)
+ update `ep_str_wrap` (2021-08-28, Sat)
+ `cnetplot` now works with a named list (2021-08-23, Mon; clusterProfiler#362)

# enrichplot 1.13.1

+ use `aplot::plot_list` instead of `cowplot::plot_grid` (2021-06-13, Sun
+ add `color_category` and `color_gene` parameters for `cnetplot`(2021-6-11, Fri)
+ Enables `showCategory` parameter to support character input in `dotplot.compareClusterResult`(2021-6-10, Thu)

# enrichplot 1.12.0

+ Bioconductor 3.13 release

# enrichplot 1.11.3

+ add function `ssplot` for similarity space plot. (2021-4-22, Thu).
+ Reconstruct the `emapplot` function and replace `emapplot_cluster` by `emapplot(group_category = TRUE)` 
+ fix bug in `emapplot_cluster.enrichResult` when the number of cluster is 2 (2021-2-24, Wed).
+ fix bug in `treeplot`: The legend is not the right size (2021-2-6, Sat).
+ fix `dotplot` for `label_format` parameter doesn't work(2021-2-3, Wed).
+ fix bug in `gseaplot2`(2021-1-28, Thu)

# enrichplot 1.11.2

+ update document (2021-1-7, Thu)
+ update `dotplot`: replace `ggsymbol::geom_symbol` with `ggstar::geom_star`(2021-1-6, Wed)
+ add parameter `shadowtext` for three functions: `emapplot`, `emapplot_cluster` and `cnetplot`. (2021-1-5, Tue)
+ update `dotplot`: supports the use of shapes and line colors to distinguish groups (2021-1-3, Sun)
+ add `treeplot` function (2020-12-29, Tue)
+ rename function `get_ww` to `get_similarity_matrix` (2020-12-29, Tue)
+ move the `emapplot` related functions to emapplot_utilities.R
+ fix bug in `emapplot` and `cnetplot` when enrichment result is one line (2020-12-26, Sat)
+ fix `pairwise_termsim` for the bug of repeated filtering of `showCategory`(2020-12-23, Wed)
+ fix `showCategory` for `cnetplot`, `emapplot`, `emapplot_cluster` when  `showCategory` is a vector of term descriptions


# enrichplot 1.11.1

+ add `orderBy` and `decreasing` parameters for `ridgeplot()` (2020-11-19, Thu)
  - <https://github.com/YuLab-SMU/enrichplot/pull/84/>
+ update `emapplot_cluster()` to label cluster in center by default and use `ggrepel` if setting `repel = TRUE` (2020-11-08, Mon)
  - <https://github.com/YuLab-SMU/enrichplot/pull/81>
+ add a `label_format` parameter to support formatting label (2020-10-28, Wed)
  + if provided with a numeric value will simply string wrap by default
  + if provided with a function will instead set labels = user_defined_function() within the scale function
  + <https://github.com/YuLab-SMU/enrichplot/pull/73>

# enrichplot 1.10.0

+ Bioconductor 3.12 release (2020-10-28, Wed)

# enrichplot 1.9.5

+ fix `wordcloud_i` (2020-10-15, Thu)
+ Remove similarity calculation from emapplot

# enrichplot 1.9.4

+ implement `pairwise_termsim` to calculate similarity of enriched terms (2020-10-09, Fri)
  - <https://github.com/YuLab-SMU/enrichplot/pull/67>
+ change parameters to be more consistent
  - <https://github.com/YuLab-SMU/enrichplot/pull/62>

# enrichplot 1.9.3

+ add `node_label_size` parameter to adjust the size of node label in `emapplot` function (2020-09-18, Fri)

# enrichplot 1.9.2

+ add function `emapplot_cluster` (2020-09-01, Tue)


# enrichplot 1.7.3

+ update `barplot` to remove using `coord_flip()` (2020-09-10, Thu)
+ update `cnetplot` color scale to tolerate with skewed foldchange (2020-03-13, Fri)
  - <https://github.com/YuLab-SMU/enrichplot/pull/40>

# enrichplot 1.7.1

+ `cnetplot` for `compareClusterResult` (`compareCluster` output) (2019-12-02, Mon)
+ move `barplot`, `dotplot` and `fortify` methods of `compareClusterResult` from `clusterProfiler` (2019-11-2, Sat)

# enrichplot 1.6.0

+ Bioconductor 3.10 release

# enrichplot 1.5.2

+ update `node_label` parameter in `cnetplot` to support selection of subset to be labeled (2019-09-27, Fri)
  - <https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#fig:cnetNodeLabel>
+ `upsetplot` for `gseaResult` (2019-09-25, Wed)
+ reimplement `upsetplot` based on `ggupset` 

# enrichplot 1.5.1

+ `gseadist` for plotting logFC distribution of selected gene sets. (2019-06-25, Tue)

# enrichplot 1.4.0

+ Bioconductor 3.9 release

# enrichplot 1.3.2

+ `dotplot` supports setting `x` to other variable, e.g. NES (2019-01-10, Thu)
+ mv vignette to [clusterProfiler-book](https://yulab-smu.github.io/clusterProfiler-book/).

# enrichplot 1.2.0

+ Bioconductor 3.8 release

# enrichplot 1.1.5

+ `gsearank` for plotting ranked list of genes belong to specific gene set
  (2018-07-04, Wed)

# enrichplot 1.1.4

+ `base_size` parameter in `gseaplot2` (2018-06-21, Thu)

# enrichplot 1.1.3

+ `pmcplot` for plotting pubmed trend (2018-06-14, Thu)
+ `ggtable` for plotting table
+ `gseaplot2` now accepts a vector of `geneSetID` (2018-06-13, Wed)

# enrichplot 1.1.2

+ `emapplot` supports `showCategory` parameter to accept a vector of
`Description`  (2018-05-29, Tue)
+ bug fixed of `showCategory` parameter for vector of `Description` in
  `cnetplot`
  - <https://support.bioconductor.org/p/109438/#109451>
+ `gseaplot2` that mimic the figure generated by broad institute's GSEA software
  (2018-05-28, Mon)

# enrichplot 1.1.1

+ `cnetplot` supports `showCategory` parameter to accept a vector of
`Description`
  - <https://github.com/GuangchuangYu/DOSE/issues/20#issuecomment-391802809>

# enrichplot 1.0.0

+ Bioconductor 3.7 release

# enrichplot 0.99.14

+ `node_label = TRUE` parameter in `cnetplot` (2018-04-08, Sun
)
+ drop NA in `dotplot` <2018-03-19, Mon>
  - <https://twitter.com/S_Canchi/status/974440351162294272>
+ enable using formula to specify x axis in `dotplot`

# enrichplot 0.99.13

+ fixed `goplot` issue by imporint `ggraph` <2018-03-12, Mon>
  - <https://github.com/GuangchuangYu/enrichplot/issues/5>

  - >Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  :
  >invalid line type
+ `dotplot` now supports `orderBy` and `decreasing` parameters to specify the order of dots by `order(x[[orderBy]], decreasing=decreasing)`


# enrichplot 0.99.9

+ defined `upsetplot` (2018-01-30, Tue)
+ all visualization methods were defined as `S4` methods (2018-01-29, Mon)

# enrichplot 0.99.5

+ defined all visualization functions as generic functions (2018-01-03, Wed)
+ add `colorEdge` parameter in `cnetplot`
+ update docs

enrichplot 0.99.3
------------------------
 + import `ggplot2::rel` to fix R check (2017-11-28, Tue)

enrichplot 0.99.0
------------------------
 + ready to submit to Bioconductor (2017-11-28, Tue)

enrichplot 0.0.3
------------------------
 + `heatplot` and `gseaplot` (2017-11-28, Tue)
 + `ridgeplot`, `barplot` and `dotplot` derived from `DOSE` (2017-11-28, Tue)
 + `cnetplot` (2017-11-28, Tue)

enrichplot 0.0.2
------------------------
 + vignette added (2017-11-28, Tue)
 + `goplot` for plotting induced GO DAG (2017-11-27, Mon)

enrichplot 0.0.1
------------------------
 + `emapplot` for plotting enrichment map (2017-11-23)
