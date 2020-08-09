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
#' @param n.core Number of cores for parallel execution (only available from Linux or Mac command line)
#' 
#' @return Data frame with maximum LOD and minimum deltaDIC for each iteration
#' 
#' @examples
#' \dontrun{
#'   par1 <- set_params(data = diallel_example,
#'                      trait = "tuber_shape")
#'                       
#'   ans1_permut <- scan1_permute(data = diallel_example,
#'                                chrom = 10,
#'                                trait = "tuber_shape",
#'                                params = par1,
#'                                n.permute = 100)
#'                              
#' }
#' 
#' @export

scan1_permute <- function(data,trait,params,n.permute=1000,chrom=NULL,dominance=1,cofactor=NULL,n.core=1) {
  LOD <- deltaDIC <- numeric(n.permute)
  for (i in 1:n.permute) {
    if (i%%100==0) {print(paste("Permutation",i))}
    data@pheno[,trait] <- sample(data@pheno[,trait])
    ans <- scan1(data=data,trait=trait,params=params,chrom=chrom,dominance=dominance,cofactor=cofactor,n.core=n.core)
    LOD[i] <- max(ans$LOD,na.rm=T)
    deltaDIC[i] <- min(ans$deltaDIC,na.rm=T)
  }
  return(data.frame(LOD=LOD,deltaDIC=deltaDIC))
}
