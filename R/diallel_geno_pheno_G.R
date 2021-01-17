#' S4 class with genotype and phenotype data and the additive relationship matrix
#' 
#' @slot ploidy Either 2 or 4
#' @slot polyorigin matrix of character strings from the genotype input file
#' @slot Xa list of matrices with the expected haplotype dosage (rows) for each parental origin genotype (columns)
#' @slot dominance Maximum dosage stored in slot \code{geno}. Integer 1-4 indicating 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. 
#' @slot X.GCA Incidence matrix for GCA effects
#' @slot map data frame with marker,chrom, position (cM and/or bp) and bin 
#' @slot geno list of length equal to the number of marker bins. Each element is a list of length ploidy corresponding to additive, digenic, trigenic, and quadrigenic effects. The elements in the nested list are sparse matrices with dimensions (id x effects). 
#' @slot pheno data frame of phenotypes
#' @slot X incidence matrix for fixed effects
#' @slot Z incidence matrix for individuals
#' @slot G1 additive relationship matrix computed by \code{\link{IBDmat}}
#' 
#' @export
diallel_geno_pheno_G <- setClass("diallel_geno_pheno_G",slots=c(G1="matrix"),contains="diallel_geno_pheno")
