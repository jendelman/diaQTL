#' @importFrom BGLR BGLR
#' @importFrom stats cor
#' 
qtl1 <- function(y,X,Z,params,geno=NULL,Xcof=NULL,X.GCA=NULL) {
  model <- "BGLR(y=y,verbose=F,response_type=params$response,burnIn=params$burnIn,nIter=params$nIter,thin=1,ETA=list(x=list(X=X,model='FIXED')"
  
  if (!is.null(X.GCA)) {
    X.GCA <- Z %*% X.GCA
    model <- paste(model,"GCA=list(X=X.GCA,model='BRR',saveEffects=FALSE)",sep=",")
  }
  
  if (!is.null(geno)) {
    g <- length(geno)
    X1 <- Z %*% geno[[1]]
    model <- paste(model,"a=list(X=X1,model='BayesC',saveEffects=TRUE)",sep=",")
    if (g>1) {
      X2 <- Z %*% geno[[2]]
      model <- paste(model,"d=list(X=X2,model='BayesC',saveEffects=TRUE)",sep=",")
      if (g>2) {
        X3 <- Z %*% geno[[3]]
        model <- paste(model,"t=list(X=X3,model='BayesC',saveEffects=TRUE)",sep=",")
        if (g>3) {
          X4 <- Z %*% geno[[4]]
          model <- paste(model,"q=list(X=X4,model='BayesC',saveEffects=TRUE)",sep=",")
        }
      }
    }
  }
  if (!is.null(Xcof)) {
    model <- paste(model,"cof=list(X=Xcof,model='BayesC')",sep=",")
  }
  model <- paste(model,"))",sep="")
  
  setwd("tmp")
  ans <- eval(parse(text=model))
  setwd("..")

  ix <- which(!is.na(y))
  if (params$response=="ordinal") {
    y2 <- as.integer(y)-1
    resid <- y2 - ans$prob[,2]
    
    yhat <- ifelse(ans$prob[ix,2] > 0.5,1,0)
    if (sd(yhat)>0) {
      R2 <- cor(yhat,y2[ix])^2*100
    } else {
      R2 <- 0
    }
  } else {
    meany <- mean(y[ix])
    R2 <- sum((ans$yHat[ix]-meany)^2)/sum((y[ix]-meany)^2)*100
    resid <- y-ans$yHat
  }
  
  result <- list(LL=ans$fit$logLikAtPostMean,DIC=ans$fit$DIC,R2=R2,resid=resid)
  return(result)
}
