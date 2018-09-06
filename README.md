<!-- README.md is generated from README.Rmd. Please edit that file -->
Visualization of Functional Enrichment Result
=============================================

[![releaseVersion](https://img.shields.io/badge/release%20version-1.0.0-green.svg?style=flat)](https://bioconductor.org/packages/enrichplot) [![develVersion](https://img.shields.io/badge/devel%20version-1.1.0-green.svg?style=flat)](https://github.com/guangchuangyu/enrichplot) [![Bioc](http://www.bioconductor.org/shields/years-in-bioc/enrichplot.svg)](https://www.bioconductor.org/packages/devel/bioc/html/enrichplot.html#since)

[![download](http://www.bioconductor.org/shields/downloads/enrichplot.svg)](https://bioconductor.org/packages/stats/bioc/enrichplot) [![total](https://img.shields.io/badge/downloads-533/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/enrichplot) [![month](https://img.shields.io/badge/downloads-152/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/enrichplot)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![platform](http://www.bioconductor.org/shields/availability/devel/treeio.svg)](https://www.bioconductor.org/packages/devel/bioc/html/treeio.html#archives) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/treeio.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/treeio/) [![Last-changedate](https://img.shields.io/badge/last%20change-2018--05--02-green.svg)](https://github.com/GuangchuangYu/treeio/commits/master)

Base classes and functions for parsing and exporting phylogenetic trees. 'treeio' supports parsing analysis findings from commonly used software packages, allows linking external data to phylogeny and merging tree data obtained from different sources. It also supports exporting phylogenetic tree with heterogeneous associated data to a single tree file.

[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=enrichplot)

### Vignettes

-   [Visualization of Functional Enrichment Result](http://bioconductor.org/packages/devel/bioc/vignettes/enrichplot/inst/doc/enrichplot.html)

Authors
-------

Guangchuang YU <https://guangchuangyu.github.io>

School of Public Health, The University of Hong Kong

[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/GuangchuangYu) [![](https://img.shields.io/badge/follow%20me%20on-微信-green.svg?style=flat)](https://guangchuangyu.github.io/blog_images/biobabble.jpg) [![](https://img.shields.io/badge/打赏-支付宝/微信-green.svg?style=flat)](https://guangchuangyu.github.io/blog_images/pay_qrcode.png)

Installation
------------

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
devtools::install_github("GuangchuangYu/enrichplot")
```

Contributing
------------

We welcome any contributions! By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CONDUCT.md).
