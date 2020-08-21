#' Frequency of haplotype dosages
#' 
#' Plot the frequency of individuals with haplotype dosage above a threshold
#' 
#' Useful for visualizing selection in selfed populations. 
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param haplotypes Names of haplotypes
#' @param dosage Dosage threshold
#' @param position Either "cM" or "bp" for plotting
#' @param chrom Names of chromosomes (default is all)
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

haplo_freq <- function(data,haplotypes,dosage,position,chrom=NULL) {
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(position %in% colnames(data@map))
  stopifnot(haplotypes %in% attr(data@geno,"haplotypes"))
  if (is.null(chrom)) {
    chrom <- unique(data@map$chrom)
  }
  nchr <- length(chrom)
  map <- data@map[data@map$chrom %in% chrom,]
  m <- nrow(map)
  
  n <- nrow(data@X.GCA)
  ans <- lapply(as.list(map$marker),function(k){haplo_get(data=data,marker=k)})
  n.haplo <- length(haplotypes)
  out <- NULL
  for (i in 1:n.haplo) {
    out <- rbind(out,data.frame(marker=map$marker,haplo=haplotypes[i],freq=sapply(ans,function(x){sum(x[,haplotypes[i]]>dosage)})/(data@ploidy*n),stringsAsFactors = F))
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
    plotme <- data.frame(x=rep(x,n.haplo),
                         y=out$freq,
                         haplo=out$haplo)
    p <- ggplot(data=plotme,aes(x=x,y=y)) + facet_wrap(~haplo,ncol=1) +
      geom_line(color="#440154") +
      ylab("Frequency") +
      theme_bw() +
      theme(text = element_text(size=13),panel.grid = element_blank()) +
      xlab(x.label) + ylim(0,NA) + geom_hline(yintercept=0,linetype=2,col="gray30")
  } else {
    col <- ifelse(as.integer(factor(map$chrom))%%2==1,"1","0")
    x <- get_x(map=data.frame(chrom=map$chrom,position=x,stringsAsFactors = F))
    plotme <- data.frame(x=rep(x,n.haplo),y=out$freq,col=rep(col,n.haplo),haplo=out$haplo)
    breaks <- (tapply(x,map$chrom,max) + tapply(x,map$chrom,min))/2
    p <- ggplot(data=plotme,aes(x=x,y=y,colour=col)) + facet_wrap(~haplo,ncol=1) +
      ylab("Frequency") +
      theme_bw() +
      scale_x_continuous(name="Chromosome",breaks=breaks,labels=chrom) +
      scale_colour_manual(values=c("#21908c","#440154"))+
      geom_hline(yintercept=0,linetype=2,col="gray30") + ylim(0,NA) + 
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
