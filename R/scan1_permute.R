#' Permutation test for scan1
#' 
#' Permutation test for scan1
#' 
#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param trait Name of trait
#' @param params List containing burnIn and nIter
#' @param n.permute Number of permutations
#' @param chrom Names of chromosomes to scan (default is all)
#' @param dominance Dominance degree (1-4)
#' @param cofactor Optional name of marker to include as cofactor in the scan
#' @param n.core Number of cores for parallel execution
#' 
#' @return Data frame with maximum LOD and minimum deltaDIC for each iteration
#' 
#' @examples
#' \dontrun{
#'   set_params(data = diallel_example,
#'              trait = "tuber_shape")
#'                       
#'   ans1_permut <- scan1_permute(data = diallel_example,
#'                                chrom = 10,
#'                                trait = "tuber_shape",
#'                                params = list(burnIn=60,nIter=600),
#'                                n.permute = 100)
#'                                
#'   ## computing permutation threshold for alpha=0.05                            
#'   quantile(ans1_permut$min_deltaDIC, 0.05)
#'                              
#' }
#' 
#' @export
#' @importFrom Matrix Matrix

scan1_permute <- function(data,trait,params,n.permute=1000,chrom=NULL,dominance=1,cofactor=NULL,n.core=1) {
  LOD <- deltaDIC <- numeric(n.permute)
  n <- nrow(data@pheno)
  for (i in 1:n.permute) {
    print(paste("Permutation",i))
    ix <- sample(1:n)
    data@pheno <- data@pheno[ix,]
    data@X <- Matrix(data@X[ix,])
    ans <- scan1(data=data,trait=trait,params=params,chrom=chrom,dominance=dominance,cofactor=cofactor,n.core=n.core)
    deltaDIC[i] <- min(ans$deltaDIC,na.rm=T)
  }
  return(data.frame(permutation=1:n.permute,min_deltaDIC=deltaDIC))
}
