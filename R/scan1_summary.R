#' Summary of scan1 result
#' 
#' Summary of scan1 result
#' 
#' @param scan1_data output from scan1
#' @param thresh optional, threshold for plotting
#' @param chrom optional, subset of chromosomes to plot 
#' @param position Either "cM" (default) or "bp"
#' @param flip should QTL be plotted as peaks (TRUE) or valleys (FALSE)
#' 
#' @return List containing 
#' \describe{
#' \item{peaks}{data frame of markers with the lowest DIC on each chromosome}
#' \item{plot}{ggplot object}
#' }
#' @examples
#' \dontrun{
#'   scan1_summary( scan1_example )
#'   scan1_summary( scan1_example, chrom = "10" )
#'   scan1_summary( scan1_example, chrom = c( "10", "12" ) ) 
#'   }
#' @export
#' @import ggplot2
#' @importFrom rlang .data

scan1_summary <- function(scan1_data,
                          thresh=NULL,
                          chrom = NULL,
                          position = "cM",
                          flip = TRUE) {

  statistic = "deltaDIC"
  stopifnot(position %in% colnames(scan1_data))
  stopifnot(statistic %in% c("LOD","deltaDIC"))
  
  if(!is.null(chrom)) {
    scan1_data <- scan1_data[scan1_data$chrom %in% chrom,]
  }
  if(position=="bp") {
    x <- scan1_data[,position]/1e6
    x.label <- "Position (Mb)"
  } else {
    x <- scan1_data[,position]
    x.label <- "Position (cM)"
  }
  
  allchr <- unique(scan1_data$chrom)
  nchr <- length(allchr)
  
  k <- integer(nchr)
  for (i in 1:nchr) {
    y <- scan1_data[,statistic]
    y[scan1_data$chrom!=allchr[i]] <- NA
    if(statistic == "deltaDIC"){
      k[i] <- which.min(y)
    }else{
      k[i] <- which.max(y)
    }
  }
  
  #if(flip){
  #  scan1_data[,"deltaDIC"] = -scan1_data[,"deltaDIC"]
    #scan1_data[,"LOD"] = -scan1_data[,"LOD"]
  #}
  # if(statistic=="deltaDIC"){
  #   statisticName = ifelse(flip, "-\U0394 DIC","\U0394 DIC")
  # }else{
  #   statisticName = ifelse(flip, "-\U0394 LOD","LOD")
  # }
  statisticName = "\U0394 DIC"
  
  p <- NULL
  if (nchr==1) {
    plotme <- data.frame(x=x,
                         y=scan1_data[,statistic],
                         chrom=scan1_data$chrom) 
    p <- ggplot(data=plotme,aes(x=.data$x,y=.data$y)) +
        geom_line(color="#440154") +
        ylab(statisticName) + scale_y_reverse() + 
        theme_bw() +
        theme(text = element_text(size=13),panel.grid = element_blank()) +
        xlab(x.label)
  } else {
    col <- ifelse(as.integer(factor(scan1_data$chrom))%%2==1,"1","0")
    x <- get_x(map=data.frame(chrom=scan1_data$chrom,position=x,stringsAsFactors = F))
    plotme <- data.frame(x=x,y=scan1_data[,statistic],col=col)
    breaks <- (tapply(x,scan1_data$chrom,max) + tapply(x,scan1_data$chrom,min))/2
    p <- ggplot(data=plotme,aes(x=.data$x,y=.data$y,colour=.data$col)) +
        ylab(statisticName) +
        theme_bw() +
        scale_x_continuous(name="Chromosome",breaks=breaks,labels=allchr) +
        scale_colour_manual(values=c("#21908c","#440154")) + scale_y_reverse() + 
        theme(text = element_text(size=13),panel.grid = element_blank(),legend.position = "none")
    
    for (i in allchr){
      ix <- which(scan1_data$chrom==i)
      p <- p + geom_line(data=plotme[ix,])
    }
  }
    
  if (!is.null(thresh)) {
    if(flip)
      thresh=-thresh
    p <- p + geom_hline(yintercept = thresh,linetype="dashed",colour="#B89600")
  }

  #peaks <- scan1_data[k,]
  #peaks$LOD <- round(peaks$LOD,1)
  #if(flip){
  #  peaks$LOD=-peaks$LOD
  #  peaks$deltaDIC=-peaks$deltaDIC
  #}
  return(list(peaks=scan1_data[k,],plot=p))
}
