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
