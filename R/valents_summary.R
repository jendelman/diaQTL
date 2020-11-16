#' Summarize valents results
#' 
#' Tabulate the proportion of bivalents and quadrivalents from PolyOrigin
#' 
#'  The valents input file should be generated with \code{\link{read_polyancestry}}.
#' 
#' @param filename Name of valents input file
#' @return List containing
#' \describe{
#' \item{result}{data frame with results}
#' \item{plot}{stacked bar chart from ggplot2}
#' }
#'  
#' @export
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#' @import ggplot2

valents_summary <- function(filename) {
  data <- read.csv(filename,as.is=T,check.names=F)
  parents <- sort(unique(c(data$parent1,data$parent2)))
  data$parent1 <- factor(data$parent1,levels=parents)
  data$parent2 <- factor(data$parent2,levels=parents)
  n.par <- length(parents)
  chroms <- unique(data$chrom)
  n.chrom <- length(chroms)
  val <- strsplit(colnames(data)[-(1:4)],split="|",fixed=T)
  par1 <- sapply(val,function(z){z[1]})
  par2 <- sapply(val,function(z){z[2]})
  for (k in 5:8) {
    par2 <- gsub(as.character(k),replacement=as.character(k-4),x=par2)
  }
  n <- nrow(data)
  par1.prob <- matrix(0,nrow=n,ncol=4)
  colnames(par1.prob) <- par1[1:4]
  par2.prob <- par1.prob
  for (i in 1:n) {
    par1.prob[i,] <- tapply(as.numeric(data[i,-(1:4)]),factor(par1,levels=par1[1:4]),sum)*data$n[i]
    par2.prob[i,] <- tapply(as.numeric(data[i,-(1:4)]),factor(par2,levels=par1[1:4]),sum)*data$n[i]
  }
  
  tmp1 <- apply(par1.prob,2,function(z){tapply(z,list(data$parent1,data$chrom),sum)})
  tmp1[is.na(tmp1)] <- 0
  tmp2 <- apply(par2.prob,2,function(z){tapply(z,list(data$parent2,data$chrom),sum)})
  tmp2[is.na(tmp2)] <- 0
  numerator <- tmp1 + tmp2
  tmp1 <- tapply(data$n,list(data$parent1,data$chrom),sum)
  tmp1[is.na(tmp1)] <- 0
  tmp2 <- tapply(data$n,list(data$parent2,data$chrom),sum)
  tmp2[is.na(tmp2)] <- 0
  denominator <- matrix(rep(as.vector(tmp1 + tmp2),4),ncol=4,byrow=F)
  
  result <- data.frame(expand.grid(parent=parents,chrom=chroms),numerator/denominator)
  tmp <- colnames(result)
  tmp[-(1:2)] <- par1[1:4]
  colnames(result) <- tmp
  plotdata <- pivot_longer(data=result,cols=3:6,names_to="valents")
  plotdata$chrom <- factor(plotdata$chrom)
  plotdata$valents <- factor(plotdata$valents,levels=par1[1:4],ordered=T)
  p <- ggplot(data=plotdata,aes(x=parent,y=value,fill=valents)) + geom_col() + facet_wrap(~chrom) + scale_fill_brewer(name="Valents",palette="Set1") + theme(axis.text.x=element_text(angle=90,vjust=0.5)) + ylab("Proportion") + xlab("")
  
  result$parent <- as.character(result$parent)
  result$chrom <- as.character(result$chrom)
  return(list(result=result,plot=p))
}
