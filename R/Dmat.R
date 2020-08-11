#' Dominance matrix
#' 
#' Calculates the dominance (D) relationship matrix from founder genotype probabilities
#' 
#' Parameter \code{dominance} refers to 2 = digenic, 3 = trigenic, 4 = quadrigenic (Gallais 2003).  Can specify to use only a subset of the chromosomes (by default, all chromosomes are used). Calculated based on the marker bins.
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param chrom Optional, vector of chromosome names to include
#' @param dominance Either 2, 3, or 4 
#' 
#' @return Dominance relationship matrix
#' 
#' @references Gallais, A. 2003. Quantitative Genetics and Breeding Methods in  Autopolyploid Plants. Institut National de la Recherche Agronomique, Paris. 
#' 
#' @examples
#' \dontrun{
#'   Dmat_example = Dmat(data = diallel_example, dominance=2) #digenic dominance
#'   Dmat_example = Dmat(data = diallel_example, dominance=3) #trigenic dominance
#' }
#' 
#' 
#' @export

Dmat <- function(data,chrom=NULL,dominance=2) {
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(dominance %in% 2:4)
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
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
  if (dominance==2) {
    coeff <- ifelse(data@ploidy==2,1,6)
  }
  if (dominance==3) {
    coeff <- 4
  }
  if (dominance==4) {
    coeff <- 1
  }
    
  for (i in map$marker) {
    K <- K + as.matrix(tcrossprod(data@geno[[i]][[dominance]]))/(coeff^2)
  }
  return(coeff*K/m)
}