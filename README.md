<!-- README.md is generated from README.Rmd. Please edit that file -->

# Visualization of Functional Enrichment Result

[![](https://img.shields.io/badge/release%20version-1.6.0-green.svg)](https://www.bioconductor.org/packages/enrichplot)
[![](https://img.shields.io/badge/devel%20version-1.7.1-green.svg)](https://github.com/guangchuangyu/enrichplot)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/enrichplot.svg)](https://www.bioconductor.org/packages/devel/bioc/html/enrichplot.html#since)

[![download](http://www.bioconductor.org/shields/downloads/release/enrichplot.svg)](https://bioconductor.org/packages/stats/bioc/enrichplot)
[![](https://img.shields.io/badge/download-69141/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/enrichplot)
[![](https://img.shields.io/badge/download-5170/month-blue.svg)](https://bioconductor.org/packages/stats/bioc/enrichplot)

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![platform](http://www.bioconductor.org/shields/availability/devel/treeio.svg)](https://www.bioconductor.org/packages/devel/bioc/html/treeio.html#archives)
[![Build
Status](http://www.bioconductor.org/shields/build/devel/bioc/treeio.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/treeio/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--12--02-green.svg)](https://github.com/GuangchuangYu/treeio/commits/master)

The ‘enrichplot’ package implements several visualization methods for
interpreting functional enrichment results obtained from ORA or GSEA
analysis. All the visualization methods are developed based on ‘ggplot2’
graphics.

For details, please visit
<https://yulab-smu.top/clusterProfiler-book/> or
<https://yulab-smu.top/biomedical-knowledge-mining-book/>.

## :writing\_hand: Authors

Guangchuang YU <https://guangchuangyu.github.io>

School of Basic Medical Sciences, Southern Medical University

[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=enrichplot)
[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/GuangchuangYu)
[![](https://img.shields.io/badge/follow%20me%20on-WeChat-green.svg)](https://guangchuangyu.github.io/blog_images/biobabble.jpg)

## :arrow\_double\_down: Installation

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
## install.packages("devtools")
devtools::install_github("YuLab-SMU/enrichplot")
```

## :sparkling\_heart: Contributing

We welcome any contributions\! By participating in this project you
agree to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
