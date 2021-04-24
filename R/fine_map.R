#' Visualize haplotype switches for fine mapping
#' 
#' Visualize haplotype switches for fine mapping
#' 
#' Function returns graphic for all individuals with a haplotype switch (defined as change in dosage from 0 to \eqn{\geq} 1 or vice versa) for \code{haplotype} within \code{interval}. If \code{trait} is included, the trait values for each individual are displayed on the right side. The function requires map positions in bp to be included in \code{data}.
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param haplotype Name of parental haplotype
#' @param interval 2-vector with marker names
#' @param trait Name of trait to plot (optional)
#' 
#' @return ggplot2 variable
#' 
#' @import ggplot2
#' @export

fine_map <- function(data,haplotype,interval,trait=NULL) {
  
  stopifnot(haplotype %in% attr(data@geno,"haplotypes"))
  stopifnot(all(interval %in% data@map$marker))
  stopifnot("bp" %in% colnames(data@map))
  if (!is.null(trait)) {
    stopifnot(trait %in% colnames(data@pheno))
  }
  
  id <- attr(data@geno,"id")
  n <- length(id)
  
  markers <- c(get_bin(interval[1],data@map),get_bin(interval[2],data@map))
  chrom <- data@map$chrom[match(markers[1],data@map$marker)]
  map <- data@map[data@map$chrom==chrom,]
  bins <- map$bin[match(markers,map$marker)]
  map <- map[map$bin %in% bins[1]:bins[2],]
  m <- nrow(map)
  
  hapans <- matrix(0,nrow=n,ncol=m)
  i <- 0
  id.ans <- integer(0)
  for (j in 1:n) {
    y <- haplo_get(data=data,id=id[j])
    y <- y[map$marker,haplotype]
    if (min(y)==0 & max(y) >= 1) {
      i <- i + 1
      hapans[i,] <- y
      id.ans <- append(id.ans,j)
    }
  }
  n2 <- length(id.ans)
  if (n2==0) {
    print("No haplotype switches in this interval")
    return(NULL)
  }
  hapans <- hapans[1:n2,]
  
  if (n2 > 1) {
    decreasing <- which(hapans[,1]>=0.5)
    if (length(decreasing) > 1) {
      tmp <- apply(hapans[decreasing,],1,function(z){min(which(z==0))})
      decreasing <- decreasing[order(tmp,decreasing = F)]
    } 
    increasing <- setdiff(1:n2,decreasing)
    if (length(increasing) > 1) {
      tmp <- apply(hapans[increasing,],1,function(z){max(which(z==0))})
      increasing <- increasing[order(tmp,decreasing = T)]
    }
    ix <- c(decreasing,increasing)
    hapans <- hapans[ix,]
    id.ans <- id.ans[ix]
  }
  
  xmin <- 1:m - 0.5
  xmax <- 1:m + 0.5
  ymin <- seq(n2-1,0,by=-1)
  ymax <- seq(n2,1,by=-1)
  plotme <- data.frame(xmin=rep(xmin,n2),xmax=rep(xmax,n2),ymin=rep(ymin,each=m),ymax=rep(ymax,each=m),z=as.numeric(t(hapans)))
  p <- ggplot(data=plotme) + geom_rect(aes(fill=z,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) + theme_bw() + scale_fill_distiller(name="Dosage",palette="Blues",direction=1) + scale_x_continuous(name="",breaks=1:m,labels=map$marker,sec.axis=dup_axis(labels=map$bp)) + theme(axis.text.x.top = element_text(angle=90,hjust=1,vjust=0.5)) + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + scale_y_continuous(name="",breaks=ymax-0.5,labels=id[id.ans]) + geom_vline(xintercept=1:m,colour="gray50",size=0.3)
  
  if (!is.null(trait)) {
    mean.y <- tapply(data@pheno[,trait],data@pheno$id,mean)
    trait.y <- mean.y[match(id[id.ans],names(mean.y))]
    p <- p + geom_label(data=data.frame(x=rep(m+3,n2),y=ymax-0.5,label=format(trait.y,digits=2)),mapping=aes(x=x,y=y,label=label)) 
  }
  return(p)
}
