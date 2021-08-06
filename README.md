# Introduction

`MicrobiomeProfiler` is an R-shiny package based on `clusterProfiler`,
is user-friendly tool for functional enrichment analysis in microbiome
data ( G Yu, et al., 2012).

# Installation

You can install `MicrobiomeProfiler` from github using devtools :

    if (!requireNamespace("devtools", quietly=TRUE))
        install.packages("devtools")
    devtools::install_github("YuLab-SMU/MicrobiomeProfiler")

# Quick Guides

    library(MicrobiomeProfiler)
    run_app()

# Reference

***G Yu***, LG Wang, Y Han, QY He. clusterProfiler: an R package for
comparing biological themes among gene clusters. *OMICS: A Journal of
Integrative Biology* 2012, 16(5):284-287.
<doi:%5B10.1089/omi.2011.0118>\](<http://dx.doi.org/10.1089/omi.2011.0118>)