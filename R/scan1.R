#' Single QTL scan
#' 
#' Performs a linear regression for each position in the map. 
#' 
#' Parameter \code{dominance} has possible values of 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. MCMC \code{params} can be estimated using \code{\link{set_params}}. Optional argument \code{covariate} is used to include other markers in the model during the scan, which can improve statistical power with multiple QTL. It is a data frame with three columns: marker = name of the marker, dominance = 1 to 4, and epistasis = TRUE/FALSE. Function returns deltaDIC = DIC for the QTL model relative to null model with only GCA effects for the parents, as well as LL = posterior mean of the log-likelihood, which is used by \code{\link{BayesCI}}. 

#' @param data variable of class \code{\link{diallel_geno_pheno}}
#' @param trait name of trait
#' @param params list containing burnIn and nIter
#' @param dominance maximum dominance for the scan, see Details
#' @param covariate optional, to include markers as covariates. See Details
#' @param chrom names of chromosomes to scan (default is all)
#' @param n.core number of cores for parallel execution
#' 
#' @return Data frame containing the map, LL, and deltaDIC. 
#' 
#' @examples
#' \dontrun{
#' ## getting minimum burnIn and nIter
#'   set_params(data = diallel_example,
#'              trait = "tuber_shape")
#'                       
#'   scan1_example <- scan1(data = diallel_example,
#'                          chrom = "10",
#'                          trait = "tuber_shape",
#'                          params = list(burnIn=60,nIter=600))
#' }
#' 
#' @export
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport

scan1 <- function(data,trait,params=list(burnIn=100,nIter=1000),dominance=1,covariate=NULL,chrom=NULL,n.core=1) {
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
  
  if (is.null(chrom)) {
    chrom <- unique(data@map$chrom)
  }

  y <- data@pheno[,trait]
  if (inherits(y,"factor")) {
    response <- "ordinal"
  } else {
    response <- "gaussian"
  }
  params <- list(response=response,nIter=max(params$nIter),burnIn=max(params$burnIn))
  
  #no marker model
  ans0 <- runBGLR(params=params,y=y,Xfix=data@X,Z=data@Z,Xgca=data@X.GCA,saveEffects=FALSE)
  
  #with marker
  map <- data@map[data@map$chrom %in% chrom,]
  bins <- get_bin(map$marker,map)

  f1 <- function(marker,data,dominance,covariate,trait,params){
    if (!is.null(covariate)) {
      qtl <- data.frame(marker=c(marker,covariate$marker),dominance=c(dominance,covariate$dominance))
      ix <- which(covariate$epistasis)
      n.epi <- length(ix)
      if (n.epi > 0) {
        epistasis <- data.frame(marker1=rep(marker,n.epi),marker2=covariate$marker[ix])
      } else {
        epistasis <- NULL
      }
    } else {
      qtl <- data.frame(marker=marker,dominance=dominance)
      epistasis <- NULL
    }
    ans <- fitQTL(data=data,trait=trait,qtl=qtl,epistasis=epistasis,
                  polygenic=FALSE,params=params,CI.prob=-1)
    return(ans)
  }
  
  cl <- makeCluster(n.core)
  clusterExport(cl=cl,varlist=NULL)
  ans1 <- parLapply(cl, unique(bins),f1,data=data,dominance=dominance,covariate=covariate,trait=trait,params=params)
  stopCluster(cl)
  
  fit <- data.frame(LL=as.numeric(sapply(ans1,function(x){x$LL})),deltaDIC=as.numeric(sapply(ans1,function(x){x$DIC})-ans0$DIC))
  return(data.frame(map[,1:(ncol(map)-1)],fit[match(bins,unique(bins)),],stringsAsFactors = F))
}
