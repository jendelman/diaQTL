diaQTL R Package
================
Jeffrey Endelman and Rodrigo Amadeu

This R package is for QTL analysis of diploid and autotetraploid diallel populations. Phenotypes are regressed on genotype probabilities, and the regression coefficients are random effects. R package BGLR is used for the regression analysis. The PolyOrigin software is recommended to calculate the genotype probabilities.

A [vignette is available](https://jendelman.github.io/diaQTL/diaQTL_vignette.html) to illustrate the workflow using a sample potato dataset. More detailed information is available in the [reference manual](https://jendelman.github.io/diaQTL/diaQTL_manual.pdf).

Financial support for this research comes from USDA NIFA Award No. 2019-67013-29166

To install and load the package:

``` r
install.packages("devtools")
library(devtools)
install_github("jendelman/diaQTL", build_vignettes=FALSE)
library(diaQTL)
```
