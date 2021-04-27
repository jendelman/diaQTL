get_trace <- function(qtl,epistasis,polygenic,params) {
  n.qtl <- nrow(qtl)
  qtl.trace <- vector("list",n.qtl)
  for (i in 1:n.qtl) {
    dominance <- qtl$dominance[i]
    qtl.trace[[i]] <- vector("list",dominance)
    for (j in 1:dominance) {
      ans <- readBinMat(sub("Q",paste(i,j,sep="_"),x="tmp/ETA_qtlQ_b.bin"))
      qtl.trace[[i]][[j]] <- t(apply(ans,1,function(z){z-mean(z)}))
    }
  }
  if (!is.null(epistasis)) {
    n.epi <- nrow(epistasis)
    epi.trace <- vector("list",n.epi)
    for (i in 1:n.epi) {
      ans <- readBinMat(sub("Q",i,x="tmp/ETA_aaQ_b.bin"))
      epi.trace[[i]] <- t(apply(ans,1,function(z){z-mean(z)}))
    }
  } else {
    epi.trace <- NULL
  }
  
  if (polygenic) {
    poly.trace <- scan("tmp/ETA_poly_varU.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]
  } else {
    poly.trace <- NULL
  }
  resid.trace <- scan("tmp/varE.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]
  return(list(qtl=qtl.trace,epi=epi.trace,poly=poly.trace,resid=resid.trace))
}