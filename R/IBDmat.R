#' Realized IBD relationship
#' 
#' Calculates realized relationship matrices (additive and dominance) from founder genotype probabilities
#' 
#' Parameter \code{dominance} refers to 1 = additive, 2 = digenic, 3 = trigenic, 4 = quadrigenic (Gallais 2003).  Can specify to use only a subset of the chromosomes (by default, all chromosomes are used). Calculated for each chromosome based on the marker bins and then averaged across chromosomes.
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param dominance One of 1,2,3,4 
#' @param chrom Optional, vector of chromosome names to include
#' @param n.core number of cores for parallel execution
#' 
#' @return Relationship matrix
#' 
#' @references Gallais, A. 2003. Quantitative Genetics and Breeding Methods in  Autopolyploid Plants. Institut National de la Recherche Agronomique, Paris. 
#' 
#' @examples
#' \dontrun{
#'   IBD_example = IBDmat(data = diallel_example, dominance=1) #additive
#'   IBD_example = IBDmat(data = diallel_example, dominance=2) #digenic dominance
#' }
#' 
#' 
#' @export

IBDmat <- function(data,dominance=1,chrom=NULL,n.core=1) {
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(dominance %in% 1:4)
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
  if (is.null(chrom)) {
    chrom <- unique(data@map$chrom)
  }
  stopifnot(all(chrom %in% data@map$chrom))
  
  bin.names <- names(data@geno)
  map <- data@map[data@map$marker %in% bin.names & data@map$chrom %in% chrom,]
  markers <- split(map$marker,f=map$chrom)
  m <- nrow(map)
  
  f1 <- function(ix,data,dominance) {
    id <- attr(data@geno,"id")
    n <- length(id)
    K <- matrix(0,nrow=n,ncol=n)
    dimnames(K) <- list(id,id)
    for (i in ix) {
      K <- K + as.matrix(tcrossprod(data@geno[[i]][[dominance]]))
    }
    m <- length(ix)
    return(K/m/choose(n=data@ploidy,k=dominance))
  }

  if (dominance > 1) {  
    cl <- makeCluster(n.core)
    clusterExport(cl=cl,varlist=NULL)
    ans <- parLapply(cl,markers,f1,data=data,dominance=dominance)
    stopCluster(cl)
  } else {
    ans <- data@A[chrom]
  }

  n.chrom <- length(chrom)
  ans2 <- ans[[1]]
  if (n.chrom > 1) {
    for (i in 2:n.chrom) {
      ans2 <- ans2 + ans[[i]]
    }
  }
  return(ans2/n.chrom)
}
