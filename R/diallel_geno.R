#' S4 class with genotype data
#' 
#' @slot ploidy Either 2 or 4
#' @slot polyorigin matrix of character strings from the genotype input file, one row per bin
#' @slot Xa list of matrices with the expected haplotype dosage (rows) for each parental origin genotype (columns)
#' @slot dominance Maximum dosage stored in slot \code{geno}. Integer 1-4 indicating 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. 
#' @slot X.GCA Incidence matrix for GCA effects
#' @slot map data frame with marker,chrom, position (cM and/or bp) and bin 
#' @slot geno list of length equal to the number of marker bins. Each element is a list of length \code{dominance}. The elements in the nested list are sparse matrices with dimensions (id x effects), containing the dosage for each effect.
#' @slot prob list of length equal to the number of marker bins. Each element is a sparse matrix (dimensions id x haplotypes) containing the probability that the haplotype is present in the individual.
#' 
#' @export
diallel_geno <- setClass("diallel_geno",slots=c(ploidy="integer",polyorigin="matrix",Xa="list",dominance="integer",X.GCA="Matrix",map="data.frame",geno="list",prob="list"))
