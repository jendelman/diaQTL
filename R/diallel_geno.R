#' S4 class with genotype data
#' 
#' @slot ploidy Either 2 or 4
#' @slot dominance Integer 1-4 indicating 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance.
#' @slot X.GCA Incidence matrix for GCA effects
#' @slot map data frame with marker,chrom, position (cM and/or bp) and bin 
#' @slot geno list of length equal to the number of marker bins. Each element is a list of length \code{dominance}. The elements in the nested list are sparse matrices with dimensions (id x effects). 
#' 
#' @export
diallel_geno <- setClass("diallel_geno",slots=c(ploidy="integer",dominance="integer",X.GCA="Matrix",map="data.frame",geno="list"))
