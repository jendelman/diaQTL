#' Summary of scan1 result
#' 
#' Summary of scan1 result
#' 
#' @param scan1_data output from scan1
#' @param display Logical variable whether to plot the LOD score
#' @param thresh optional, LOD threshold for plotting
#' @param chromosome string with chrom name(s) to plot. By default, all chromosomes are plotted
#' @param distance string with "cM" for centiMorgans or "bp" for basepairs
#' 
#' @return List containing 
#' \describe{
#' \item{peaks}{Data frame of the markers with the highest LOD score per chromosome}
#' \item{plot}{ggplot object}
#' }
#' @examples
#' \dontrun{
#'   scan1_summary( scan1_example )
#'   scan1_summary( scan1_example, chromosome = "10" )
#'   scan1_summary( scan1_example, chromosome = c( "10", "12" ) ) 
#'   }
#' @export
#' @import ggplot2

scan1_summary <- function(scan1_data,
                          display=T,
                          thresh=NULL,
                          chromosome = NULL,
                          distance = "cM") {
  if(!is.null(chromosome)) {
    scan1_data <- scan1_data[scan1_data$chromosome %in% chromosome,]
  }
  if(distance=="bp") {
    scan1_data$position <- scan1_data$position/1e6
    x.label <- "Position (Mb)"
  } else {
    scan1_data$position <- scan1_data$position
    x.label <- "Position (cM)"
  }
  
  allchr <- unique(scan1_data$chromosome)
  nchr <- length(allchr)
  
  k <- integer(nchr)
  for (i in 1:nchr) {
    y <- scan1_data$LOD
    y[scan1_data$chromosome!=allchr[i]] <- NA
    k[i] <- which.max(y)
  }
  
  p <- NULL
  if (display) {
    #make a plot
    if (nchr==1) {
      plotme <- data.frame(x=scan1_data$position,
                           y=scan1_data$LOD,
                           chrom=scan1_data$chromosome) 
      p <- ggplot(data=plotme,aes(x=x,y=y)) +
        geom_line(color="#440154") +
        ylab("LOD") +
        theme_bw() +
        theme(text = element_text(size=13),panel.grid = element_blank()) +
        xlab(x.label)
    } else {
      col <- ifelse(as.integer(factor(scan1_data$chromosome))%%2==1,"1","0")
      x <- get_x(scan1_data[,c("chromosome","position")])
      plotme <- data.frame(x=x,y=scan1_data$LOD,col=col)
      breaks <- (tapply(x,scan1_data$chromosome,max) + tapply(x,scan1_data$chromosome,min))/2
      p <- ggplot(data=plotme,aes(x=x,y=y,colour=col)) +
        ylab("LOD") +
        theme_bw() +
        scale_x_continuous(name="Chromosome",breaks=breaks,labels=allchr) +
        scale_colour_manual(values=c("#21908c","#440154"))+
        theme(text = element_text(size=13),panel.grid = element_blank(),legend.position = "none")
    
      for (i in allchr){
        ix <- which(scan1_data$chromosome==i)
        p <- p + geom_line(data=plotme[ix,])
      }
    }
    
    if (!is.null(thresh)) {
      p <- p + geom_hline(yintercept = thresh,linetype="dashed",colour="#B89600")
    }
  }
  
  peaks <- scan1_data[k,-match("position",colnames(scan1_data))]
  peaks$R2 <- round(peaks$R2,1)
  peaks$LOD <- round(peaks$LOD,1)
  peaks$deltaDIC <- round(peaks$deltaDIC,1)
  return(list(peaks=peaks,plot=p))
}
