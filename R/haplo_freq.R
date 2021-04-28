#' Haplotype frequencies
#' 
#' Plots the frequency of individuals with haplotype dosage above a threshold
#' 
#' Useful for visualizing selection in selfed populations. For multiple chromosomes, each haplotype is shown in its own panel using facet_wrap. For one chromosome, the haplotypes are shown on the same set of axes.
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param haplotypes Names of haplotypes
#' @param dosage Dosage threshold
#' @param id Vector of id names (default is entire population)
#' @param position Either "cM" (default) or "bp" for plotting
#' @param chrom Names of chromosomes (default is all)
#' @param markers Optional, markers to indicate with dashed line. Only available when plotting a single chromosome.
#' 
#' @return List containing
#' \describe{
#' \item{result}{Data frame with the map and frequency}
#' \item{plot}{ggplot object}
#' }
#' 
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data

haplo_freq <- function(data,haplotypes,dosage,id=NULL,position="cM",chrom=NULL,markers=NULL) {
  #y=haplo=z=ymin=ymax=NULL #to avoid NOTE while doing R check
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(position %in% colnames(data@map))
  stopifnot(haplotypes %in% attr(data@geno,"haplotypes"))
  stopifnot(all(id %in% rownames(data@X.GCA)))
  if (is.null(chrom)) {
    chrom <- unique(data@map$chrom)
  }
  nchr <- length(chrom)
  map <- data@map[data@map$chrom %in% chrom,]
  m <- nrow(map)
  
  if (is.null(id)) {
    id  <- rownames(data@X.GCA)
  }
  n <- length(id)
  ans <- lapply(as.list(map$marker),function(k){haplo_get(data=data,marker=k)})
  n.haplo <- length(haplotypes)
  out <- NULL
  for (i in 1:n.haplo) {
    out <- rbind(out,data.frame(marker=map$marker,haplo=haplotypes[i],freq=sapply(ans,function(x){sum(x[id,haplotypes[i]]>dosage)})/n,stringsAsFactors = F))
  }
  out$haplo <- factor(out$haplo)

  if(position=="bp") {
    x <- map$bp/1e6
    x.label <- "Position (Mb)"
  } else {
    x <- map$cM
    x.label <- "Position (cM)"
  }
  
  if (nchr==1) {
    features <- markers
    plotme <- data.frame(x=rep(x,n.haplo),
                         y=out$freq,
                         haplo=out$haplo)
    p <- ggplot(data=plotme,aes(x=.data$x,y=.data$y,colour=.data$haplo)) + geom_line() + scale_colour_brewer(name="Haplotype",palette="Set1") + ylab("Frequency") + theme_bw() + theme(text = element_text(size=13)) + xlab(x.label) + ylim(0,1) + geom_hline(yintercept=0,linetype=2,col="gray30") + geom_hline(yintercept=1,linetype=2,col="gray30")
    if (!is.null(features)) {
      p <- p + labs(subtitle = paste(c("Markers:",features),collapse=" "))
      for (q in 1:length(features)) {
        marker <- features[q]
        if (position=="cM") {
          marker <- get_bin(marker=marker,map=data@map)
        }
        stopifnot(marker %in% map$marker)
        k <- match(marker,map$marker)
        p <- p + geom_vline(xintercept=x[k],linetype=2)
      }
    } 
  } else {
    col <- ifelse(as.integer(factor(map$chrom))%%2==1,"1","0")
    x <- get_x(map=data.frame(chrom=map$chrom,position=x,stringsAsFactors = F))
    plotme <- data.frame(x=rep(x,n.haplo),y=out$freq,col=rep(col,n.haplo),haplo=out$haplo)
    breaks <- (tapply(x,map$chrom,max) + tapply(x,map$chrom,min))/2
    p <- ggplot(data=plotme,aes(x=.data$x,y=.data$y,colour=.data$col)) + facet_wrap(~.data$haplo,ncol=1) +
      ylab("Frequency") +
      theme_bw() +
      scale_x_continuous(name="Chromosome",breaks=breaks,labels=chrom) +
      scale_colour_manual(values=c("#21908c","#440154"))+
      geom_hline(yintercept=0,linetype=2,col="gray30") + ylim(0,1) + 
      theme(text = element_text(size=13),panel.grid = element_blank(),legend.position = "none")
    for (q in 1:n.haplo) {
      for (i in chrom){
        ix <- which(map$chrom==i)
        p <- p + geom_line(data=plotme[ix+(q-1)*m,])
      }
    }
  }
  
  tmp <- as.data.frame(pivot_wider(data=out,names_from="haplo",values_from="freq"))
  result <- merge(map[,1:(ncol(map)-1)],tmp)
  result <- result[match(map$marker,result$marker),]
  return(list(result=result,plot=p))
}
