#' Single QTL scan
#' 
#' Performs a linear regression for each position in the map. 
#' 
#' LOD score is the difference between the log10-likelihood for the QTL model vs. no QTL model (higher is better). deltaDIC is the difference between the Deviance Information Criterion for the QTL model vs. no QTL model (lower values is better). r2 is the squared correlation between the fitted and observed values. Parameter \code{dominance} controls the genetic model: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. MCMC \code{params} can be estimated using \code{\link{set_params}}.  

#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param trait Name of trait
#' @param params List containing burnIn and nIter
#' @param dominance Dominance degree (1-4)
#' @param chrom Names of chromosomes to scan (default is all)
#' @param cofactor Optional name of marker to include as cofactor in the scan
#' @param n.core Number of cores for parallel execution (only available from Linux or Mac command line)
#' 
#' @return Data frame containing the map, LOD, r2 and deltaDIC results. 
#' 
#' @examples
#' \dontrun{
#'   par1 <- set_params(data = diallel_example,
#'                      trait = "tuber_shape")
#'                       
#'   scan1_example <- scan1(data = diallel_example,
#'                 chrom = 10,
#'                 trait = "tuber_shape",
#'                 params = par1)
#' }
#' 
#' @export
#' @importFrom BGLR BGLR
#' @importFrom stats model.matrix
#' @importFrom parallel mclapply

scan1 <- function(data,trait,params,dominance=1,chrom=NULL,cofactor=NULL,n.core=1) {
  
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
  params <- list(response=response,nIter=params$nIter,burnIn=params$burnIn)
  
  if (!is.null(cofactor)) {
    stopifnot(cofactor %in% data@map$marker)
    Xcof <- data@Z%*%data@geno[[get_bin(cofactor,data@map)]][[1]]
  } else {
    Xcof <- NULL
  }

  #no marker model
  ans0 <- qtl1(y=y,X=data@X,Z=data@Z,params=params,Xcof=Xcof,X.GCA=data@X.GCA)
  
  #with marker
  map <- data@map[data@map$chrom %in% chrom,]
  bins <- get_bin(map$marker,map)
  ans <- mclapply(X=unique(bins),FUN=function(k,y,data,Xcof,params,dominance){qtl1(y=y,X=data@X,Z=data@Z,Xcof=Xcof,params=params,geno=data@geno[[k]][1:dominance])},y=y,data=data,Xcof=Xcof,params=params,dominance=dominance,mc.cores = n.core)
  
  fit <- data.frame(LOD=as.numeric(sapply(ans,function(x){x$LL})-ans0$LL)/log(10),r2=as.numeric(sapply(ans,function(x){x$r2})),deltaDIC=as.numeric(sapply(ans,function(x){x$DIC})-ans0$DIC))
  return(data.frame(map[,1:(ncol(map)-1)],fit[match(bins,unique(bins)),],stringsAsFactors = F))
}
