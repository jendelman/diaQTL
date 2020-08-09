#' S4 class with genotype and phenotype data
#' 
#' @slot ploidy Either 2 or 4
#' @slot dominance Integer 1-4 indicating 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance.
#' @slot X.GCA Incidence matrix for GCA effects
#' @slot map data frame with marker,chrom, position (cM and/or bp) and bin 
#' @slot geno list of length equal to the number of marker bins. Each element is a list of length ploidy corresponding to additive, digenic, trigenic, and quadrigenic effects. The elements in the nested list are sparse matrices with dimensions (id x effects). 
#' @slot pheno data frame of phenotypes
#' @slot X incidence matrix for fixed effects
#' @slot Z incidence matrix for individuals
#' 
#' @export
diallel_geno_pheno <- setClass("diallel_geno_pheno",slots=c(pheno="data.frame",X="Matrix",Z="Matrix"),contains="diallel_geno")
