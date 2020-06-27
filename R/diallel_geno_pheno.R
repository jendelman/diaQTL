#' S4 class with genotype and phenotype data
#' 
#' @slot ploidy Either 2 or 4
#' @slot ped data frame with pedigree information. Variables are id,population,mother,father
#' @slot map data frame with marker,chrom,and position (either bp or cM)
#' @slot geno list of length 2. The first element (named "A") is a list of sparse matrices, one for each marker, with dimensions (id x alleles), containing the allele dosages, which are the regression variables for additive effects. The second (optional) element (named "D") has the same structure (list of sparse matrices), but each matrix contains the dosage of allele-pairs, which are the regression variables for digenic dominance effects. 
#' @slot pheno data frame of phenotypes
#' @slot X incidence matrix for fixed effects
#' @slot Z incidence matrix for individuals
#' 
#' @export
diallel_geno_pheno <- setClass("diallel_geno_pheno",slots=c(pheno="data.frame",X="Matrix",Z="Matrix"),contains="diallel_geno")
