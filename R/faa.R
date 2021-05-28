#' Construct incidence matrix for additive x additive epistasis
#' 
#' Construct incidence matrix for additive x additive epistasis
#' 
#' @param j first marker
#' @param k second marker
#' @param data variable inheriting from class \code{\link{diallel_geno}}
#' 
#' @return Matrix with dimensions nind x (n.hap)^2
#' 
#' @keywords internal
#' @import Matrix
#' 
faa <- function(j,k,data) {
  # j is position of first marker
  # k is position of second marker
  
  tmp1 <- unlist(lapply(data@input[j,],strsplit,split="=>",fixed=T),recursive = F)
  tmp1 <- lapply(tmp1,strsplit,split="|",fixed=T)
  tmp2 <- unlist(lapply(data@input[k,],strsplit,split="=>",fixed=T),recursive = F)
  tmp2 <- lapply(tmp2,strsplit,split="|",fixed=T)
  
  states1 <- lapply(tmp1,function(x){as.integer(x[[1]])})
  states2 <- lapply(tmp2,function(x){as.integer(x[[1]])})
  
  genoprob1 <- lapply(tmp1,function(x){y<-as.numeric(x[[2]])
  y/sum(y)})
  genoprob2 <- lapply(tmp2,function(x){y<-as.numeric(x[[2]])
  y/sum(y)})
  
  tmp3 <- mapply(FUN=function(u,s1,s2,w1,w2){
    if (length(s1)==1) {
      u1 <- Matrix(u[,s1],ncol=1)
    } else {
      u1 <- u[,s1]
    }
    if (length(s2)==1) {
      u2 <- Matrix(u[,s2],ncol=1)
    } else {
      u2 <- u[,s2]
    }
    uu <- kronecker(u1,u2)
    ww <- Matrix(kronecker(w1,w2),ncol=1)
    return(uu%*%ww)
  },u=data@Xa,s1=states1,s2=states2,w1=genoprob1,w2=genoprob2)
  
  if(is.list(tmp3)){ 
    tmp3 <- sapply(tmp3,function(x){as.vector(x)})
  } 
  return(Matrix(t(tmp3)))
}
  