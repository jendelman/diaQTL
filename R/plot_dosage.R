#' Plot parental allele dosage
#' 
#' Plot parental allele dosages across the chromosome for one individual
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param id Name of individual
#' @param chrom Name of chromosome
#' @param position Either "cM" or "bp"
#' @param marker Optional, position of marker indicated with dashed line
#' 
#' @return ggplot object 
#' 
#' @examples
#' \dontrun{
#' plot_dosage(data = diallel_example, 
#'             id = "W15263-8R", 
#'             chrom = 10)
#'             
#' plot_dosage(data = diallel_example, 
#'             id = "W15263-8R", 
#'             chrom = 10,
#'             marker = "solcap_snp_c2_25522")
#' }
#' @export
#' @import ggplot2

plot_dosage <- function(data,id,chrom,distance,marker=NULL) {
  
  ix <- which(data@map[,2] == chrom)
  position <- data@map[ix,3] #x axis values
  m <- length(position)
  if (distance=="bp") {
    position <- position/1e6
    x.label <- "Position (Mb)"
  } else {
    x.label <- "Position (cM)"
  }
  xmin <- c(position[1],position[-m] + diff(position)/2)
  xmax <- c(xmin[-1],position[m])

  geno <- dosage(data,id=id)
  geno <- geno[ix,]
  alleles <- colnames(geno)
  tmp <- strsplit(split=".",x=alleles,fixed=T)
  founders <- sapply(tmp,function(x){x[1]})
  parents <- unique(c(as.character(data@ped[id,]$mother),as.character(data@ped[id,]$father)))
  n.par <- length(parents)
  iq <- which(founders %in% parents[1])
  plotme <- data.frame(z=as.vector(geno[,iq]),xmin=xmin,xmax=xmax,ymin=rep(0:3,each=m),ymax=rep(1:4,each=m))
  
  if (n.par==2) {
    #F1
    iq <- which(founders %in% parents[2])
    plotme <- rbind(plotme,data.frame(z=as.vector(geno[,iq]),xmin=xmin,xmax=xmax,ymin=rep(4:7,each=m),ymax=rep(5:8,each=m)))
  }
  p <- ggplot(data=plotme) + geom_rect(aes(fill=z,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) + theme_bw() + scale_fill_distiller(name="Dosage",palette="Blues",direction=1) + scale_y_continuous(name="",labels=alleles[which(founders %in% parents)],breaks=(1:(4*n.par))-0.5) + xlab(x.label)
  if (!is.null(marker)) {
    stopifnot(marker %in% data@map[ix,1])
    k <- match(marker,data@map[ix,1])
    p <- p + geom_vline(xintercept=position[k],linetype=2)
  }
  return(p)
}
