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
