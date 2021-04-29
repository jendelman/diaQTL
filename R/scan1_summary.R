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
#'   scan1_summary( scan1_example, chrom = "10", thresh = -20)
#'   }
#' @export
#' @import ggplot2
#' @importFrom rlang .data

scan1_summary <- function(scan1_data,
                          thresh=NULL,
                          chrom = NULL,
                          position = "cM",
                          flip = FALSE) {

  stopifnot(position %in% colnames(scan1_data))
  
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
    y <- scan1_data$deltaDIC
    y[scan1_data$chrom!=allchr[i]] <- NA
    k[i] <- which.min(y)
  }

  p <- NULL
  if (nchr==1) {
    plotme <- data.frame(x=x,
                         y=scan1_data$deltaDIC,
                         chrom=scan1_data$chrom) 
    p <- ggplot(data=plotme,aes(x=.data$x,y=.data$y)) +
        geom_line(color="#440154") +
        ylab("\U0394 DIC") + 
        theme_bw() +
        theme(text = element_text(size=13),panel.grid = element_blank()) +
        xlab(x.label)
  } else {
    col <- ifelse(as.integer(factor(scan1_data$chrom))%%2==1,"1","0")
    x <- get_x(map=data.frame(chrom=scan1_data$chrom,position=x,stringsAsFactors = F))
    plotme <- data.frame(x=x,y=scan1_data$deltaDIC,col=col)
    breaks <- (tapply(x,scan1_data$chrom,max) + tapply(x,scan1_data$chrom,min))/2
    p <- ggplot(data=plotme,aes(x=.data$x,y=.data$y,colour=.data$col)) +
        ylab("\U0394 DIC") +
        theme_bw() +
        scale_x_continuous(name="Chromosome",breaks=breaks,labels=allchr) +
        scale_colour_manual(values=c("#21908c","#440154")) + 
        theme(text = element_text(size=13),panel.grid = element_blank(),legend.position = "none")
    
    for (i in allchr){
      ix <- which(scan1_data$chrom==i)
      p <- p + geom_line(data=plotme[ix,])
    }
  }
  
  
  if (!is.null(thresh)) {
    p <- p + geom_hline(yintercept = thresh,linetype="dashed",colour="#B89600")
  }
  
  if(!flip){
    p <- p + scale_y_reverse()
  }
  return(list(peaks=scan1_data[k,],plot=p))
}
