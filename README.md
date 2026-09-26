<!-- README.md is generated from README.Rmd. Please edit that file -->

# Visualization of Functional Enrichment Result

[![](https://img.shields.io/badge/release%20version-1.32.0-green.svg)](https://www.bioconductor.org/packages/enrichplot)
[![](https://img.shields.io/badge/devel%20version-1.99.6-green.svg)](https://github.com/YuLab-SMU/enrichplot)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/enrichplot.svg)](https://www.bioconductor.org/packages/devel/bioc/html/enrichplot.html#since)

[![download](http://www.bioconductor.org/shields/downloads/release/enrichplot.svg)](https://bioconductor.org/packages/stats/bioc/enrichplot)
[![](https://img.shields.io/badge/download-1405847/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/enrichplot)
[![](https://img.shields.io/badge/download-28045/month-blue.svg)](https://bioconductor.org/packages/stats/bioc/enrichplot)

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![platform](http://www.bioconductor.org/shields/availability/devel/enrichplot.svg)](https://www.bioconductor.org/packages/devel/bioc/html/enrichplot.html#archives)
[![Build
Status](http://www.bioconductor.org/shields/build/devel/bioc/enrichplot.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/enrichplot/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2026--09--26-green.svg)](https://github.com/YuLab-SMU/enrichplot/commits/master)

The ‘enrichplot’ package provides visualization methods for interpreting
functional enrichment results from ORA or GSEA analyses. It is designed
to work with the ‘clusterProfiler’ ecosystem and builds on ‘ggplot2’ for
flexible and extensible graphics.

For details, please visit
<https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html>.

## :writing_hand: Authors

Guangchuang YU <https://yulab-smu.top>

School of Basic Medical Sciences, Southern Medical University

## :arrow_double_down: Installation

Get the released version from Bioconductor:

``` r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("enrichplot")
```

Or the development version from github:

``` r
## install.packages("remotes")
remotes::install_github("YuLab-SMU/enrichplot")

## or
## install.packages("yulab.utils")
yulab.utils::install_zip_gh("YuLab-SMU/enrichplot")
```

## :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
