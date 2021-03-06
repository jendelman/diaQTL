#' Cluster parental haplotypes
#' 
#' Cluster parental haplotypes
#' 
#' The argument \code{marker} can be either a single marker or vector of two markers. 
#' If a single marker, the function finds the smallest interval containing that marker such that 
#' the phased SNP haplotypes are all unique. If two markers are provided, that interval is used. 
#' Clustering utilizes hclust(method="average"). See also \code{\link{phased_parents}} for an additional visualization tool.
#' 
#' @param filename Name of diaQTL_parents input file
#' @param marker Either target marker or marker interval (see Details).
#' @param haplotypes Vector of haplotype names (default is all)
#' 
#' @return List containing 
#' \describe{
#' \item{haplo}{Data frame of haplotypes}
#' \item{dendro}{Dendrogram}
#' }
#' 
#' @export
#' @importFrom utils read.csv write.csv
#' @importFrom stats dist as.dendrogram hclust
#' @importFrom labeling extended
#' @importFrom ggfittext geom_fit_text
#' @import ggdendro

haplo_cluster <- function(filename,marker,haplotypes=NULL) {
  #y=yend=xend=NULL #to avoid NOTE while doing R check
  
  map <- read.csv(filename,as.is=T,check.names=F)
  stopifnot(marker %in% map$marker)
  stopifnot(length(marker) %in% 1:2)
  parents <- colnames(map)[-(1:4)]
  n.parents <- length(parents)
  all.hap <- character(0)
  ploidy <- length(gregexpr(pattern="|",map[1,5],fixed=T)[[1]])+1
  for (i in 1:n.parents) {
    all.hap <- c(all.hap,paste(parents[i],1:ploidy,sep="."))
  }
  if (!is.null(haplotypes)) {
    stopifnot(haplotypes %in% all.hap)
  } else {
    haplotypes <- all.hap
  }
  
  if (length(marker)==1) {
    map <- map[map$chrom==map$chrom[match(marker,map$marker)],]
    m <- nrow(map)
    k <- match(marker,map$marker)
    x <- expand.grid(first=1:k,last=k:m)
    x$size <- x$last-x$first
    x <- x[order(x$size,x$first),]
    x <- x[-1,]
  } else {
    iq <- match(marker,map$marker)
    map <- map[iq[1]:iq[2],]
    m <- nrow(map)
  }
  
  iu <- match(parents,colnames(map))
  if (n.parents==1) {
    tmp <- strsplit(map[,iu],split="|",fixed=T)
    ans <- matrix(as.integer(unlist(tmp)),nrow=m,ncol=ploidy,byrow=T)
  } else {
    tmp <- apply(map[,iu],2,strsplit,split="|",fixed=T)
    ploidy <- length(tmp[[1]][[1]])
    tmp2 <- lapply(tmp,function(z){matrix(as.integer(unlist(z)),ncol=ploidy,nrow=m,byrow=T)})
    ans <- NULL
    for (i in 1:n.parents) {
      ans <- cbind(ans,tmp2[[i]])
    }
  }
  colnames(ans) <- all.hap
  
  if (length(marker)==1) {
    finished <- FALSE
    j <- 1
    while (!finished & j <= nrow(x)) {
      first <- x$first[j]
      last <- x$last[j]
      d <- as.matrix(dist(t(ans[first:last,haplotypes]),method="manhattan"))
      diag(d) <- NA
      smallest <- min(apply(d,1,min,na.rm=T))
    
      if (smallest > 0) {
        finished <- TRUE
      } else {
        j <- j + 1
      }
    }
  } else {
    first <- 1
    last <- nrow(ans)
    finished <- TRUE
  }
  d <- dist(t(ans[first:last,haplotypes]),method="manhattan")
  dhc <- as.dendrogram(hclust(d,method="average"))
  ddat <- dendro_data(dhc,type="rectangle")
  xmax <- max(ddat$segments$yend)
  breaks <- extended(0,xmax,5)
  p <- ggplot(segment(ddat)) + geom_segment(aes(x=.data$y,y=.data$x,xend=.data$yend,yend=.data$xend)) + scale_y_continuous(name="",breaks=NULL,labels=NULL) + geom_fit_text(data=ddat$labels,mapping=aes(xmax=.data$y,xmin=-10,y=.data$x,label=.data$label),place="right") + scale_x_continuous("Distance",limits=c(-10,xmax),breaks=breaks,labels=breaks) + theme_classic() +  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

  if (!finished) {
    print("Reached end of the chromosome and failed to find unique haplotypes")
  }
  return(list(haplo=cbind(map[first:last,1:4],ans[first:last,]),dendro=p))
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       