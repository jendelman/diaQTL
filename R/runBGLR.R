#' Xqtl is list of length equal to the number of QTL, and each element is a list of maximum length 4, for additive up to quadrigenic dominance. Xaa is list of matrices for epistasis.
#'
#' @importFrom BGLR BGLR
#' @keywords internal
runBGLR <- function(params,y,Xfix,Z,Xgca=NULL,Xqtl=NULL,Xaa=NULL,polyG=NULL,saveEffects) {
  ix <- which(!is.na(y))
  model <- "BGLR(y=y[ix],verbose=F,response_type=params$response,burnIn=params$burnIn,nIter=params$nIter,thin=1,ETA=list(fix=list(X=Xfix[ix,],model='FIXED')"
  
  if (!is.null(Xgca)) {
    model <- paste(model,"GCA=list(X=as.matrix(Z[ix,] %*% Xgca),model='BRR',saveEffects=FALSE)",sep=",")
  }
  
  if (!is.null(Xqtl)) {
    n.qtl <- length(Xqtl)
    for (i in 1:n.qtl) {
      dominance <- length(Xqtl[[i]])
      for (j in 1:dominance) {
        model <- paste(model,gsub("j",j,gsub("Q",i,"qtlQ_j=list(X=as.matrix(Z[ix,] %*% Xqtl[[Q]][[j]]),model='BayesC',saveEffects=saveEffects)")),sep=",")
      }
    }
  }
  if (!is.null(Xaa)) {
    n.aa <- length(Xaa)
    for (i in 1:n.aa) {
      model <- paste(model,gsub("Q",i,"aaQ=list(X=as.matrix(Z[ix,] %*% Xaa[[Q]]),model='BayesC',saveEffects=saveEffects)"),sep=",")
    }
  }
  if (!is.null(polyG)) {
    model <- paste(model,"poly=list(K=as.matrix(Z[ix,] %*% polyG %*% t(Z[ix,])),model='RKHS')",sep=",")
  }
  
  model <- paste(model,"))",sep="")
  
  ans <- try(setwd("tmp"),silent=T)
  if (inherits(ans,"try-error")) {
    dir.create("tmp")
    setwd("tmp")
  } 
  ans <- suppressWarnings(eval(parse(text=model)))
  setwd("..")

  resid <- numeric(length(y))*NA
  if (params$response=="ordinal") {
    y2 <- as.integer(y[ix])-1
    resid[ix] <- y2 - ans$prob[,2]
    #yhat <- ifelse(ans$prob[,2] > 0.5,1,0)
  } else {
    y2 <- y[ix]
    resid[ix] <- y2 - ans$yHat
    #yhat <- ans$yHat[ix]
  }

  result <- list(LL=ans$fit$postMeanLogLik,DIC=ans$fit$DIC,resid=resid)
  return(result)
}
