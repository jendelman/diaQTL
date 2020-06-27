#' A matrix
#' 
#' Calculates the additive (A) relationship matrix from founder genotype probabilities
#' 
#' Additive relationships are calculated from kinship coefficients of order 2 (Gallais 2003). Can be subset to one or more chromosomes (which is useful for the leave-one-chromosome-out kinship method), or specific markers can be excluded after QTL detection.
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param chrom Only use markers from these chromosomes
#' @param exclude.marker Markers to exclude
#' 
#' @return A matrix
#' 
#' @references Gallais, A. 2003. Quantitative Genetics and Breeding Methods in  Autopolyploid Plants. Institut National de la Recherche Agronomique, Paris. 
#' 
#' @examples
#' \dontrun{
#'   Amat_example = Amat(data = diallel_example)
#'   Amat_example = Amat(data = diallel_example, chrom=c(1:11)) #leave chromosome 12 out
#'   Amat_example = Amat(data = diallel_example, exclude.marker = "solcap_snp_c2_25522")
#' }
#' 
#' @export
Amat <- function(data,chrom=NULL,exclude.marker=NULL) {
  stopifnot(inherits(data,"diallel_geno"))
  if (is.null(chrom)) {
    chrom <- unique(data@map[,2])
  }
  stopifnot(all(chrom %in% data@map[,2]))
  stopifnot(all(exclude.marker %in% data@map[,1]))
  
  ix <- which(data@map[,1] %in% exclude.marker | !(data@map[,2] %in% chrom)) #exclude
  m <- nrow(data@map) 
  m2 <- m-length(ix)
  
  n <- nrow(data@ped)
  Ka <- matrix(0,nrow=n,ncol=n)
  colnames(Ka) <- rownames(Ka) <- attr(data@geno$A,"id")
  
  for (i in setdiff(1:m,ix)) {
    Ka <- Ka + as.matrix(tcrossprod(data@geno$A[[i]]))/(data@ploidy^2)
  }
  return(A=data@ploidy*Ka/m2)
}
