#' Determine parameters for scan1
#' 
#' Determine parameters for scan1
#' 
#' The burn-in and total number of iterations are determined using the Raftery and Lewis diagnostic from R package \code{coda}, based on a 95% probability that the estimated median of the additive effects is between the quantiles (0.5-tol) to (0.5+tol). For greater precision, decrease the \code{tol} parameter. 
#' 
#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param trait Name of trait
#' @param tol tolerance for estimating the median
#' @param burnIn initial value for burnIn parameter
#' @param nIter initial value for nIter parameter
#' @param dominance dominance degree (1-4)
#' 
#' @return List containing 
#' \describe{
#' \item{burnIn}{Number of burn-in iterations}
#' \item{nIter}{Total number of iterations}
#' }
#' 
#' @examples
#' \dontrun{
#'   par1 <- set_params(data = diallel_example,
#'                      trait = "tuber_shape")
#' }
#'   
#' @export
#' @importFrom coda raftery.diag mcmc.list mcmc
#' @importFrom BGLR readBinMat
#' @importFrom stats quantile

set_params <- function(data,trait,tol=0.1,dominance=1,burnIn=50,nIter=1000) {
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
  
  y <- data@pheno[,trait]
  if (inherits(y,"factor")) {
    response <- "ordinal"
  } else {
    response <- "gaussian"
  }
  params <- list(response=response,nIter=nIter,burnIn=burnIn)
  
  chrom <- unique(data@map[,2])
  n.chrom <- length(chrom)
  raftans <- NULL
  for (i in 1:n.chrom) {
    k <- data@map$marker[match(chrom[i],data@map[,2])] #first marker for each chromosome
    ans <- qtl1(y=y,X=data@X,params=params,Z=data@Z,geno=data@geno[[k]][1:dominance])  
    tmp <- mcmc(readBinMat('ETA_add_b.bin'))
    raftans <- rbind(raftans,raftery.diag(tmp,q=0.5,r=tol)$resmatrix[,1:2])
  }
  
  return(list(burnIn=round(as.numeric(quantile(raftans[,1],probs=0.9)),0),nIter=round(as.numeric(quantile(raftans[,2],probs=0.9)),0)))
}
