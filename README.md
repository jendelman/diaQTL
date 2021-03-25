diaQTL R Package
================
Jeffrey Endelman and Rodrigo Amadeu

This R package is for QTL and haplotype analysis of autotetraploid (and
diploid) diallel populations. Phenotypes are regressed on genotype
probabilities, and the regression coefficients are random effects. R
package BGLR is used for the regression analysis. The [PolyOrigin
software](https://github.com/chaozhi/PolyOrigin.jl) is recommended to
calculate the genotype probabilities for tetraploids.

A [vignette is
available](https://jendelman.github.io/diaQTL/diaQTL_vignette.html) to
illustrate the workflow using a sample potato dataset, and more detailed
information is available in the [reference
manual](https://jendelman.github.io/diaQTL/diaQTL_manual.pdf).

Please cite our
[manuscript](https://www.biorxiv.org/content/10.1101/2020.12.18.423479v1)
if you use the package. Financial support for this research comes from
USDA NIFA Award No. 2019-67013-29166

To install and load the package:

``` r
install.packages("devtools")
library(devtools)
install_github("jendelman/diaQTL", build_vignettes=FALSE)
library(diaQTL)
```

Video resources:

1) Tutorial about `diaQTL` presented as part of the "Tools for Polyploids Workshop 2021":

[![Watch the video](https://img.youtube.com/vi/iOxckvAWCnU/maxresdefault.jpg)](https://youtu.be/iOxckvAWCnU)

2) Seminar on "QTL Mapping in Tetraploid Diallel Populations" presented as part of the Computational Genetics Discussion Group (The Roslin Institute) Seminar Series of 2021:

[![Watch the video](https://img.youtube.com/vi/UF8UJkGl16Y/maxresdefault.jpg)](https://youtu.be/UF8UJkGl16Y)
