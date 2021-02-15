#' Plot parental haplotype dosage
#' 
#' Plot parental haplotype dosages across the chromosome for one individual
#' 
#' For "cM" plotting, only one marker per bin is displayed. For "bp" plotting, all markers are included. If multiple individuals are included in \code{id}, then the plot shows the probability that a haplotype is present in all individuals. 
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param id Name of individual(s)
#' @param chrom Name of chromosome
#' @param position Either "cM" (default) or "bp"
#' @param markers Optional, markers to indicate with dashed line
#' 
#' @return ggplot object 
#' 
#' @examples
#' \dontrun{
#' haplo_plot(data = diallel_example, 
#'             id = "W15263-8R", 
#'             chrom = 10)
#'             
#' haplo_plot(data = diallel_example, 
#'             id = "W15263-8R", 
#'             chrom = 10,
#'             marker = "solcap_snp_c2_25522")
#' }
#' @export
#' @import ggplot2

haplo_plot <- function(data,id,chrom,position="cM",markers=NULL) {
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(position %in% colnames(data@map))
  stopifnot(chrom %in% data@map$chrom)
  features <- markers
  
  marks <- data@map$marker[data@map$chrom %in% chrom]
  if (position=="bp") {
    x <- data@map[data@map$marker %in% marks,"bp"] #x axis values
    x <- x/1e6
    x.label <- "Position (Mb)"
  } else {
    marks <- intersect(names(data@geno),marks) #only plot the bins
    x <- data@map[data@map$marker %in% marks,"cM"] #x axis values
    x.label <- "Position (cM)"
  }
  m <- length(x)
  xmin <- c(0,x[-m] + diff(x)/2)
  xmax <- c(xmin[-1],x[m])
  
  haplotypes <- attr(data@geno,"haplotypes")
  tmp <- strsplit(split=".",x=haplotypes,fixed=T)
  founders <- sapply(tmp,function(x){x[1]})
  parents <- lapply(id,function(id){unique(colnames(data@X.GCA)[data@X.GCA[id,] > 0])})
  stopifnot(sapply(parents,function(z){all(z==parents[[1]])}))  #parents must be the same
  
  parents <- parents[[1]]
  n.par <- length(parents)
  iq <- which(founders %in% parents[1])
  y1 <- 0:(data@ploidy-1)
  
  geno <- haplo_get(data,id=id[1])
  n <- length(id)
  if (n > 1) {
    for (k in 2:n) {
      geno <- geno*haplo_get(data,id=id[k])
    }
  }
  geno <- geno[marks,]
  plotme <- data.frame(z=as.vector(geno[,iq]),xmin=xmin,xmax=xmax,ymin=rep(y1,each=m),ymax=rep(y1+1,each=m))
  
  if (n.par==2) {
    #F1
    y2 <- y1+data@ploidy
    iq <- which(founders %in% parents[2])
    plotme <- rbind(plotme,data.frame(z=as.vector(geno[,iq]),xmin=xmin,xmax=xmax,ymin=rep(y2,each=m),ymax=rep(y2+1,each=m)))
  }
  
  p <- ggplot(data=plotme) + geom_rect(aes(fill=z,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) + theme_bw() + scale_fill_distiller(name="Dosage",palette="Blues",direction=1) + scale_y_continuous(name="",labels=haplotypes[which(founders %in% parents)],breaks=(1:(data@ploidy*n.par))-0.5) + xlab(x.label) + geom_hline(yintercept=0:(data@ploidy*n.par),color="gray30")
  id <- paste(id,collapse=" ")
  if (!is.null(features)) {
    p <- p + labs(title = id,subtitle = paste(c("Markers:",features),collapse=" "))
    for (q in 1:length(features)) {
      marker <- features[q]
      if (position=="cM") {
        marker <- get_bin(marker=marker,map=data@map)
      }
      stopifnot(marker %in% marks)
      k <- match(marker,marks)
      p <- p + geom_vline(xintercept=x[k],linetype=2)
    }
  } else {
    p <- p + labs(title=id)
  }
  return(p)
}
