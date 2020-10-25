#' Determine burn-in and total number of iterations
#' 
#' Determine burn-in and total number of iterations
#' 
#' Determines the burn-in and total number of iterations using the Raftery and Lewis diagnostic from R package \code{coda}, based on a 95% probability that the estimate for quantile \code{q} of the additive effects is within the interval \code{(q-r,q+r)}. The 90th percentile for burn-in and total iterations across the additive effects is returned. Parameter \code{dominance} specifies which genetic model (1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance) to use when determining the number of iterations, but this must be parameter must still be specified when calling functions such as \code{\link{scan1}} or \code{\link{fitQTL}}. Suggested values for \code{\link{scan1}} are q=0.5 and r=0.1. For \code{\link{fitQTL}}, the values depend on the desired Bayesican credible interval. For a 90% CI, suggested values are q=0.05 and r=0.05. If \code{marker=NULL} (default), the first marker of every chromosome is analyzed to generate parameters suitable for \code{\link{scan1}}. Parameter \code{nIter} sets the number of iterations used to apply the Raftery and Lewis diagnostic; the default value is 2000, and if a larger number is needed, an error will be generated with this information.
#' 
#' @param data variable of class \code{\link{diallel_geno_pheno}}
#' @param trait name of trait
#' @param marker name of marker (optional)
#' @param dominance dominance degree 
#' @param q quantile to estimate
#' @param r tolerance for quantile
#' @param nIter number of iterations 
#' 
#' @return List containing 
#' \describe{
#' \item{burnIn}{Number of burn-in iterations}
#' \item{nIter}{Total number of iterations}
#' }
#' 
#' @examples
#' \dontrun{
#'   # Parameters for scan1
#'   par1 <- set_params(data = diallel_example,
#'                      trait = "tuber_shape",q=0.5,r=0.1)
#'                      
#'   # Parameters for fitQTL
#'   par2 <- set_params(data = diallel_example,
#'                      trait = "tuber_shape",q=0.05,r=0.05,marker="solcap_c2_25522")
#' }
#'   
#' @export
#' @importFrom coda raftery.diag mcmc.list mcmc
#' @importFrom BGLR readBinMat
#' @importFrom stats quantile

set_params <- function(data,trait,dominance=1,marker=NULL,q=0.5,r=0.1,nIter=2000) {
  
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
  params <- list(response=response,nIter=nIter,burnIn=0)  
  
  raftans <- NULL
  if (is.null(marker)) {
    chrom <- unique(data@map[,2])
    marker <- data@map$marker[match(chrom,data@map[,2])] #first marker of each chromosome
  } 
  marker <- get_bin(marker,data@map)
  n.mark <- length(marker)
  
  for (i in 1:n.mark) {
    k <- marker[i]
    ans <- qtl1(y=y,X=data@X,params=params,Z=data@Z,geno=data@geno[[k]][1:dominance])  
    tmp <- mcmc(readBinMat("tmp/ETA_a1_b.bin"))
    tmp2 <- raftery.diag(tmp,q=q,r=r)$resmatrix
    if (inherits(tmp2,"character")) {
      stop(paste("Try increasing the number of iterations to",as.integer(tmp2[2])+100))
    }
    raftans <- rbind(raftans,tmp2[,1:2])
  }
  
  return(list(burnIn=round(as.numeric(quantile(raftans[,1],probs=0.9)),0),nIter=round(as.numeric(quantile(raftans[,2],probs=0.9)),0)))
}
