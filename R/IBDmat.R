#' Realized IBD relationship
#' 
#' Calculates realized relationship matrices from founder genotype probabilities
#' 
#' Parameter \code{dominance} refers to 1 = additive, 2 = digenic, 3 = trigenic, 4 = quadrigenic.  Can specify to use only a subset of the chromosomes (by default, all chromosomes are used). If \code{epistasis} is TRUE, then \code{dominance} must be 1 (additive x additive epistasis). Only pairs of markers on different chromosomes are used for epistasis. 
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param dominance One of 1,2,3,4 
#' @param epistasis TRUE/FALSE
#' @param spacing spacing between marker bins, in cM
#' @param chrom Optional, vector of chromosome names to include
#' @param n.core number of cores for parallel execution
#' 
#' @return Relationship matrix
#' 
#' @examples
#' \dontrun{
#'   IBD_example = IBDmat(data = diallel_example, dominance=1) #additive
#'   IBD_example = IBDmat(data = diallel_example, dominance=2) #digenic dominance
#'   IBD_example = IBDmat(data = diallel_example, epistasis=TRUE) #additive x additive epistasis
#' }
#' 
#' 
#' @export
#' @importFrom parallel makeCluster clusterExport stopCluster parLapply

IBDmat <- function(data,dominance=1,epistasis=FALSE,
                   spacing=1,chrom=NULL,n.core=1) {
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(dominance %in% 0:4)
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
  if ((dominance > 1) & epistasis) {
    stop("Only additive x additive epistasis is available")
  }
  if (is.null(chrom)) {
    chrom <- unique(data@map$chrom)
  }
  stopifnot(all(chrom %in% data@map$chrom))
  
  bin.names <- names(data@geno)
  map <- data@map[data@map$marker %in% bin.names & data@map$chrom %in% chrom,]
  markers <- split(map$marker,f=map$chrom)
  pos <- split(map$cM,f=map$chrom)
  max.cM <- ceiling(tapply(map$cM,map$chrom,max))
  tmp <- mapply(FUN=function(x,max.cM){cut(x,breaks=seq(0,max.cM,by=spacing),include.lowest=T)},
                x=pos,max.cM=max.cM)
  if (!is.list(tmp)) {
    tmp <- as.list(data.frame(tmp))
  }
  ix <- lapply(tmp,function(x){match(unique(x),x)})
  markers <- mapply(FUN=function(x,ix){x[ix]},x=markers,ix=ix) 
  if (!is.list(markers)) {
    markers <- as.list(data.frame(markers))
  }
  names(markers) <- names(pos)
  map <- map[map$marker %in% unlist(markers),]
  
  if (epistasis) {
    tmp <- expand.grid(marker1=map$marker,marker2=map$marker)
    tmp$marker1 <- as.character(tmp$marker1)
    tmp$marker2 <- as.character(tmp$marker2)
    tmp$chr1 <- map$chrom[match(tmp$marker1,map$marker)]
    tmp$chr2 <- map$chrom[match(tmp$marker2,map$marker)]
    tmp <- tmp[tmp$chr1 > tmp$chr2,]
    marker.pairs <- split(tmp[,1:2],apply(tmp[,3:4],1,paste,collapse=" "))
  } 
  
  f1 <- function(x,data,dominance) {
    id <- attr(data@geno,"id")
    n <- length(id)
    K <- matrix(0,nrow=n,ncol=n)
    dimnames(K) <- list(id,id)
    m <- length(x)
    for (i in x) {
      K <- K + as.matrix(tcrossprod(data@geno[[i]][[dominance]]))
    }
    return(K/m/choose(n=data@ploidy,k=dominance))
  }

  f2 <- function(x,data) {
    id <- attr(data@geno,"id")
    n <- length(id)
    K <- matrix(0,nrow=n,ncol=n)
    dimnames(K) <- list(id,id)
    m <- nrow(x)
    for (i in 1:m) {
      K <- K + as.matrix(tcrossprod(faa(j=x$marker1[i],k=x$marker2[i],data=data)))
    }
    return(K/m/data@ploidy^2)
  }
    
  if (dominance > 0) {
    cl <- makeCluster(n.core)
    clusterExport(cl=cl,varlist=NULL)
    if (epistasis) {
      ans <- parLapply(cl,marker.pairs,f2,data=data)
    } else {
      ans <- parLapply(cl,markers,f1,data=data,dominance=dominance)
    }
    stopCluster(cl)
  } else {
    #undocumented option, used by fitQTL
    ans <- data@A[chrom]
  }
  
  k <- length(ans)
  ans2 <- ans[[1]]
  if (k > 1) {
    for (i in 2:k) {
      ans2 <- ans2 + ans[[i]]
    }
  }
  return(ans2/k)
}
