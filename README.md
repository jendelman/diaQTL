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
available](https://jendelman.github.io/diaQTL/diaQTL_Vignette.html) to
illustrate the workflow using a sample potato dataset, and more detailed
information is available in the [reference
manual](https://jendelman.github.io/diaQTL/diaQTL_manual.pdf).

Please [cite our manuscript](https://doi.org/10.1093/genetics/iyab124) if you use the package. Financial support for this 
research comes from USDA NIFA Award No. 2019-67013-29166

To install and load the package:

``` r
install.packages("devtools")
devtools::install_github("jendelman/diaQTL", build_vignettes=FALSE)
library(diaQTL)
```
