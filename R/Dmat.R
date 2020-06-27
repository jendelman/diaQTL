#' D matrix
#' 
#' Calculates the dominance (D) relationship matrix from founder genotype probabilities
#' 
#' Dominance relationships are calculated from kinship coefficients of order 4 (Gallais 2003). An entire chromosome or specific markers can be excluded.
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param chrom Only use markers from these chromosomes
#' @param exclude.marker Markers to exclude
#' 
#' @return Digenic dominance relationship matrix
#' 
#' @references Gallais, A. 2003. Quantitative Genetics and Breeding Methods in  Autopolyploid Plants. Institut National de la Recherche Agronomique, Paris. 
#' 
#' @examples
#' \dontrun{
#'   Dmat_example = Dmat(data = diallel_example)
#'   Dmat_example = Dmat(data = diallel_example, chrom=c(10,11)) #leave chromosome 12 out
#'   Dmat_example = Dmat(data = diallel_example, exclude.marker = "solcap_snp_c2_25522")
#' }
#' 
#' 
#' @export

Dmat <- function(data,chrom=NULL,exclude.marker=NULL) {
  stopifnot(inherits(data,"diallel_geno"))
  if (is.null(data@geno$D)) {
    stop("Dominance was FALSE in read_data")
  }
  if (is.null(chrom)) {
    chrom <- unique(data@map[,2])
  }
  stopifnot(all(chrom %in% data@map[,2]))
  stopifnot(all(exclude.marker %in% data@map[,1]))
  
  ix <- which(data@map[,1] %in% exclude.marker | !(data@map[,2] %in% chrom)) #exclude
  m <- nrow(data@map)
  m2 <- m-length(ix)
  
  n <- nrow(data@ped)
  Kd <- matrix(0,nrow=n,ncol=n)
  colnames(Kd) <- rownames(Kd) <- attr(data@geno$D,"id")
  coeff <- ifelse(data@ploidy==2,1,6)
    
  for (i in setdiff(1:m,ix)) {
    Kd <- Kd + tcrossprod(data@geno$D[[i]])/(coeff^2)
  }
  return(D=coeff*Kd/m2)
}