#' Single QTL scan
#' 
#' Performs a linear regression for each position in the map. 
#' 
#' For non-binary traits, R2 is the proportion of variance explained by the regression. 
#' For binary traits, R2 is the squared phi correlation.
#' LOD score is the difference 
#' between the log-likelihood of the model with QTL and the null model
#' (no QTL); higher values are better. deltaDIC is the difference between the DIC of the model with QTL minus
#' the DIC of the null model; lower values are better. 

#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param trait Name of trait
#' @param params List containing burnIn and nIter, use set_params function to estimate it
#' @param chrom Names of chromosomes to scan (default is all)
#' @param dominance Logical variable whether to include digenic dominance effects
#' @param cofactor Optional name of marker to include as cofactor in the scan
#' @param n.core Number of cores to use for parallel execution by forking (do not use with GUI)
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
#' @importFrom parallel mcmapply

scan1 <- function(data,trait,params,chrom=NULL,dominance=F,cofactor=NULL,n.core=1) {
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  stopifnot(cofactor %in% names(data@geno$A))
  if (is.null(data@geno$D) & dominance) {
    stop("Dominance was FALSE in read_data")
  }

  if (is.null(chrom)) {
    chrom <- unique(data@map[,2])
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
    Xcof <- data@Z%*%data@geno$A[[cofactor]]
  } else {
    Xcof <- NULL
  }

  #no marker model
  ans0 <- qtl1(y=y,X=data@X,Z=data@Z,params=params,Xcof=Xcof)
  
  fit <- NULL
  for (i in 1:n.chrom) {
    cat(paste("Chromosome",chrom[i],"\n"))
    ix <- which(data@map[,2]==chrom[i])
    m <- length(ix)
    
    if (dominance) {
      ans1 <- mcmapply(FUN=qtl1,genoA=data@geno$A[ix],genoD=data@geno$D[ix],MoreArgs=list(y=y,X=data@X,Z=data@Z,Xcof=Xcof,params=params),SIMPLIFY=FALSE,mc.cores=n.core)
    } else {
      ans1 <- mcmapply(FUN=qtl1,genoA=data@geno$A[ix],MoreArgs=list(y=y,X=data@X,Z=data@Z,Xcof=Xcof,params=params),SIMPLIFY=FALSE,mc.cores=n.core)
    }
    fit <- rbind(fit,data.frame(data@map[ix,],LOD=as.numeric(sapply(ans1,function(x){x$LL})-ans0$LL)/log(10),R2=as.numeric(sapply(ans1,function(x){x$R2})),deltaDIC=as.numeric(sapply(ans1,function(x){x$DIC})-ans0$DIC),stringsAsFactors = F))
  }
  return(fit)
}
