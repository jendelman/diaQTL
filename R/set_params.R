#' Determine number of iterations for MCMC
#' 
#' Determine number of iterations for MCMC
#' 
#' Determines the burn-in and total number of iterations using the Raftery and Lewis diagnostic from R package \code{coda}, based on a 95% probability that the estimate for quantile \code{q} of the additive genetic variance is within the interval \code{(q-r,q+r)}. If \code{marker=NULL} (default), the first marker of each chromosome is analyzed, and the largest value across this set is returned. Parameter \code{dominance} specifies which genetic model (1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance) to use when determining the number of iterations, but this parameter must still be specified when calling functions such as \code{\link{scan1}} or \code{\link{fitQTL}}. The default values of q=0.5 and r=0.1 are recommended for \code{\link{scan1}} based on the idea of estimating the posterior mean. For estimating the 90% Bayesian CI with \code{\link{fitQTL}}, suggested values are q=0.05, r=0.025. Parameter \code{nIter} sets the number of iterations used to apply the Raftery and Lewis diagnostic; the default value is 2000, and if a larger number is needed, an error will be generated with this information. 
#' 
#' @param data variable of class \code{\link{diallel_geno_pheno}}
#' @param trait name of trait
#' @param qtl optional data frame, see Details
#' @param epistasis optional data frame, see Details
#' @param polygenic TRUE/FALSE whether to include additive polygenic effect
#' @param q quantile to estimate
#' @param r tolerance for quantile
#' @param nIter number of iterations 
#' 
#' @return matrix showing the number of burn-in and total iterations for the genetic variances in the model
#' 
#' @examples
#' \dontrun{
#'   # Parameters for scan1
#'   par1 <- set_params(data = diallel_example,
#'                      trait = "tuber_shape",q=0.5,r=0.1)
#'                      
#'   # Parameters for fitQTL
#'   par2 <- set_params(data = diallel_example,
#'                      trait = "tuber_shape",q=0.05,r=0.025,marker="solcap_c2_25522")
#' }
#'   
#' @export
#' @importFrom coda raftery.diag mcmc
#' @importFrom BGLR readBinMat

set_params <- function(data,trait,qtl=NULL,epistasis=NULL,polygenic=FALSE,q=0.5,r=0.1,nIter=2000){
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  
  if (is.null(qtl)) {
    chrom <- unique(data@map$chrom)
    markers <- data@map$marker[match(chrom,data@map$chrom)] #first marker of each chromosome
    n.fit <- length(markers)
  } else {
    n.fit <- 1
    stopifnot(qtl$marker %in% data@map$marker)
    if (max(qtl$dominance) > data@dominance) {
      stop("Dominance degree cannot exceed value used with read_data.")
    }
  }
  
  y <- data@pheno[,trait]
  if (inherits(y,"factor")) {
    response <- "ordinal"
  } else {
    response <- "gaussian"
  }
  params <- list(response=response,nIter=nIter,burnIn=-1)  

  raftans <- NULL
  for (i in 1:n.fit) {
    if (is.null(qtl)) {
      qtl <- data.frame(marker=markers[i],dominance=1)
    }
    ans <- fitQTL(data=data,trait=trait,qtl=qtl,epistasis=epistasis,polygenic=polygenic,
                  params=params,CI.prob=NULL)
    tmp2 <- raftery.diag(mcmc(ans),q=q,r=r)$resmatrix
    if (inherits(tmp2,"character")) {
      stop(paste("Try increasing the number of iterations to",as.integer(tmp2[2])+100))
    }
    tmp2 <- matrix(tmp2[,1:2],ncol=2)
    colnames(tmp2) <- c("burnIn","nIter")
    rownames(tmp2) <- colnames(ans)
    if (is.null(raftans) || max(raftans[,2]) < max(tmp2[,2])) {
      raftans <- tmp2
    }
  }
  return(raftans)
}
