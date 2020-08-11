#' A matrix
#' 
#' Calculates the additive (A) relationship matrix from founder genotype probabilities
#' 
#' Additive relationships are calculated from kinship coefficients of order 2 (Gallais 2003). Can specify to use only a subset of the chromosomes (by default, all chromosomes are used). Calculated based on the marker bins.
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param chrom Optional, vector of chromosome names to include
#' 
#' @return A matrix
#' 
#' @references Gallais, A. 2003. Quantitative Genetics and Breeding Methods in  Autopolyploid Plants. Institut National de la Recherche Agronomique, Paris. 
#' 
#' @examples
#' \dontrun{
#'   Amat_example = Amat(data = diallel_example)
#'   Amat_example = Amat(data = diallel_example, chrom=c(1:11)) #leave chromosome 12 out
#' }
#' 
#' @export
Amat <- function(data,chrom=NULL) {
  stopifnot(inherits(data,"diallel_geno"))
  if (is.null(chrom)) {
    chrom <- unique(data@map$chrom)
  }
  stopifnot(all(chrom %in% data@map$chrom))
  bin.names <- names(data@geno)
  map <- data@map[data@map$marker %in% bin.names & data@map$chrom %in% chrom,]
  m <- nrow(map)
  
  id <- attr(data@geno,"id")
  n <- length(id)
  K <- matrix(0,nrow=n,ncol=n)
  colnames(K) <- rownames(K) <- id
  for (i in map$marker) {
    K <- K + as.matrix(tcrossprod(data@geno[[i]][[1]]))/(data@ploidy^2)
  }
  return(data@ploidy*K/m)
}
