#' @importFrom BGLR BGLR
#' @importFrom stats cor
#' 
qtl1 <- function(y,X,Z,params,geno=NULL,Xcof=NULL,G1=NULL,X.GCA=NULL) {
  model <- "BGLR(y=y,verbose=F,response_type=params$response,burnIn=params$burnIn,nIter=params$nIter,thin=1,ETA=list(x=list(X=X,model='FIXED')"
  
  if (!is.null(X.GCA)) {
    X.GCA <- Z %*% X.GCA
    model <- paste(model,"GCA=list(X=X.GCA,model='BRR',saveEffects=FALSE)",sep=",")
  }
  
  if (!is.null(geno)) {
    g <- length(geno)
    for (i in 1:g) {
      eval(parse(text=gsub("Q",i,"XQ <- Z %*% geno[[Q]]")))
      model <- paste(model,gsub("Q",i,"aQ=list(X=XQ,model='BayesC',saveEffects=TRUE)"),sep=",")
    }
  }
  if (!is.null(G1)) {
    model <- paste(model,"polyg=list(K=G1,model='RKHS')",sep=",")
  }
  if (!is.null(Xcof)) {
    model <- paste(model,"cof=list(X=Xcof,model='BayesC')",sep=",")
  }
  model <- paste(model,"))",sep="")
  
  ans <- try(setwd("tmp"),silent=T)
  if (inherits(ans,"try-error")) {
    dir.create("tmp")
    setwd("tmp")
  } 
  ans <- suppressWarnings(eval(parse(text=model)))
  setwd("..")

  ix <- which(!is.na(y))
  if (params$response=="ordinal") {
    y2 <- as.integer(y)-1
    resid <- y2 - ans$prob[,2]
    yhat <- ifelse(ans$prob[ix,2] > 0.5,1,0)
  } else {
    y2 <- y
    resid <- y2 - ans$yHat
    yhat <- ans$yHat[ix]
  }
  if (sd(yhat)>0) {
    r2 <- cor(yhat,y2[ix])^2
  } else {
    r2 <- 0
  }
  
  result <- list(LL=ans$fit$logLikAtPostMean,DIC=ans$fit$DIC,r2=r2,resid=resid)
  return(result)
}
