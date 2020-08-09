#' Single QTL scan
#' 
#' Performs a linear regression for each position in the map. 
#' 
#' For non-binary traits, R2 is the proportion of variance explained by the regression. For binary traits, R2 is the squared phi correlation.
#' LOD score is the difference between the log-likelihood of the model with QTL and the null model
#' (no QTL); higher values are better. deltaDIC is the difference between the DIC of the model with QTL minus
#' the DIC of the null model; lower values are better. Parameter \code{dominance} controls the genetic model: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. Parameter \code{params} can be estimated using \code{\link{set_params}}.

#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param trait Name of trait
#' @param params List containing burnIn and nIter
#' @param dominance Dominance degree (1-4). See Details.
#' @param chrom Names of chromosomes to scan (default is all)
#' @param cofactor Optional name of marker to include as cofactor in the scan
#' @param n.core Number of cores for parallel execution (only available from Linux or Mac command line)
#' 
#' @return Data frame containing the map, LOD, R2 and deltaDIC results. 
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
  n.chrom <- length(chrom)

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
  
  bin.names <- names(data@geno)
  map <- data@map[(data@map$chrom %in% chrom) & (data@map$marker %in% bin.names),1:(ncol(data@map)-1)]
  
  #with marker
  ans <- mclapply(X=map$marker,FUN=function(k,y,data,Xcof,params,dominance){qtl1(y=y,X=data@X,Z=data@Z,Xcof=Xcof,params=params,geno=data@geno[[k]][1:dominance])},y=y,data=data,Xcof=Xcof,params=params,dominance=dominance,mc.cores = n.core)
  
  fit <- data.frame(map,LOD=as.numeric(sapply(ans,function(x){x$LL})-ans0$LL)/log(10),R2=as.numeric(sapply(ans,function(x){x$R2})),deltaDIC=as.numeric(sapply(ans,function(x){x$DIC})-ans0$DIC),stringsAsFactors = F)
  return(fit)
}
