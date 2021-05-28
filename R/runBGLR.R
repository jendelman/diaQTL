#' Xqtl is list of length equal to the number of QTL, and each element is a list of maximum length 4, for additive up to quadrigenic dominance. Xaa is list of matrices for epistasis.
#'
#' @importFrom BGLR BGLR
#' @keywords internal
runBGLR <- function(params,y,Xfix,Z,Xgca=NULL,Xqtl=NULL,Xaa=NULL,polyG=NULL,saveEffects) {
  model <- "BGLR(y=y,verbose=F,response_type=params$response,burnIn=params$burnIn,nIter=params$nIter,thin=1,ETA=list(fix=list(X=Xfix,model='FIXED')"
  
  if (!is.null(Xgca)) {
    model <- paste(model,"GCA=list(X=Z %*% Xgca,model='BRR',saveEffects=FALSE)",sep=",")
  }
  
  if (!is.null(Xqtl)) {
    n.qtl <- length(Xqtl)
    for (i in 1:n.qtl) {
      dominance <- length(Xqtl[[i]])
      for (j in 1:dominance) {
        model <- paste(model,gsub("j",j,gsub("Q",i,"qtlQ_j=list(X=Z %*% Xqtl[[Q]][[j]],model='BayesC',saveEffects=saveEffects)")),sep=",")
      }
    }
  }
  if (!is.null(Xaa)) {
    n.aa <- length(Xaa)
    for (i in 1:n.aa) {
      model <- paste(model,gsub("Q",i,"aaQ=list(X=Z %*% Xaa[[Q]],model='BayesC',saveEffects=saveEffects)"),sep=",")
    }
  }
  if (!is.null(polyG)) {
    model <- paste(model,"poly=list(K=Z %*% polyG %*% t(Z),model='RKHS')",sep=",")
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

  result <- list(LL=ans$fit$postMeanLogLik,DIC=ans$fit$DIC,resid=resid)
  return(result)
}
