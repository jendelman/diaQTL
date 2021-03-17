#' Single QTL scan
#' 
#' Performs a linear regression for each position in the map. 
#' 
#' LOD score is the difference between the log10-likelihood for the QTL model vs. no QTL model (higher is better). deltaDIC is the difference between the Deviance Information Criterion for the QTL model vs. no QTL model (lower values is better). r2 is the squared correlation between the fitted and observed values. Parameter \code{dominance} controls the genetic model: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. MCMC \code{params} can be estimated using \code{\link{set_params}}. The argument \code{cofactor} should be a list with three components: marker = name of the marker; dominance = 1, 2, 3, or 4; epistasis = TRUE/FALSE. When a cofactor is included, the LOD and deltaDIC values are relative to a model with the cofactor. 

#' @param data variable of class \code{\link{diallel_geno_pheno}}
#' @param trait name of trait
#' @param params list containing burnIn and nIter
#' @param dominance dominance degree (1-4)
#' @param chrom names of chromosomes to scan (default is all)
#' @param cofactor optional, see Details for format.
#' @param n.core number of cores for parallel execution
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
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport

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
    stopifnot(is.list(cofactor))
    stopifnot(cofactor$marker %in% data@map$marker)
    cofactor$marker <- get_bin(cofactor$marker,data@map)
    cofactor$X <- data@geno[[cofactor$marker]][1:cofactor$dominance]
    cofactor$Xaa <- NULL
  } 

  #no marker model (but with cofactor)
  ans0 <- qtl1(y=y,X=data@X,Z=data@Z,params=params,cofactor=cofactor,X.GCA=data@X.GCA)
  
  #with marker
  map <- data@map[data@map$chrom %in% chrom,]
  bins <- get_bin(map$marker,map)

  cl <- makeCluster(n.core)
  clusterExport(cl=cl,varlist=NULL)
  parqtl1 <- function(k,y,data,params,dominance,cofactor){
    if (!is.null(cofactor) && cofactor$epistasis) {
      cofactor$Xaa <- data@Z%*%faa(j=cofactor$marker,k=k,data=data)
    }
    ans <- qtl1(y=y,X=data@X,Z=data@Z,params=params,geno=data@geno[[k]][1:dominance],cofactor)
    return(ans)
  }
  ans <- parLapply(cl, unique(bins), parqtl1, y=y, data=data, params=params, dominance=dominance, cofactor=cofactor)
  #ans <- lapply(unique(bins), parqtl1, y=y, data=data, cofactor=cofactor, params=params, dominance=dominance)
  stopCluster(cl)
  
  fit <- data.frame(LOD=as.numeric(sapply(ans,function(x){x$LL})-ans0$LL)/log(10),r2=as.numeric(sapply(ans,function(x){x$r2})),deltaDIC=as.numeric(sapply(ans,function(x){x$DIC})-ans0$DIC))
  return(data.frame(map[,1:(ncol(map)-1)],fit[match(bins,unique(bins)),],stringsAsFactors = F))
}
