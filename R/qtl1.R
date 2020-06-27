#' @importFrom BGLR BGLR
#' @importFrom stats cor
#' 
qtl1 <- function(y,X,Z,params,genoA=NULL,genoD=NULL,Xcof=NULL) {
  model <- "BGLR(y=y,verbose=F,response_type=params$response,burnIn=params$burnIn,nIter=params$nIter,thin=1,ETA=list(x=list(X=X,model='BRR')"
  
  if (!is.null(genoA)) {
    Xadd <- Z %*% genoA
    model <- paste(model,"add=list(X=Xadd,model='BayesC',saveEffects=TRUE)",sep=",")
  }
  if (!is.null(genoD)) {
    Xdom <- Z %*% genoD
    model <- paste(model,"dom=list(X=Xdom,model='BayesC',saveEffects=TRUE)",sep=",")
  }
  if (!is.null(Xcof)) {
    model <- paste(model,"cof=list(X=Xcof,model='BayesC')",sep=",")
  }
  model <- paste(model,"))",sep="")
  
  ans <- eval(parse(text=model))

  ix <- which(!is.na(y))
  if (params$response=="ordinal") {
    y2 <- as.integer(y)-1
    resid <- y2 - ans$prob[,2]
    
    yhat <- ifelse(ans$prob[ix,2] > 0.5,1,0)
    if (sd(yhat)>0) {
      R2 <- cor(yhat,y2[ix])^2*100
    } else {
      R2 <- NA
    }
  } else {
    meany <- mean(y[ix])
    R2 <- sum((ans$yHat[ix]-meany)^2)/sum((y[ix]-meany)^2)*100
    resid <- y-ans$yHat
  }
  
  result <- list(LL=ans$fit$logLikAtPostMean,DIC=ans$fit$DIC,R2=R2,resid=resid)
  return(result)
}
