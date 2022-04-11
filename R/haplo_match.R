#' Match up parental haplotypes
#' 
#' Match up parental haplotypes
#' 
#' Designed to match up parental haplotypes between two phased parent files based on genetic distance. In the plots, the haplotypes in file1 are numbered 1-4, and those in file2 are numbered 5-8.
#' 
#' @param file1 name of CSV phased parent file 
#' @param file2 name of CSV phased parent file
#' @param chrom chromosome
#' 
#' @return Data frame with results
#' @importFrom stats cutree hclust dist
#' @export
#' 
haplo_match <- function(file1, file2, chrom) {
  data1 <- read.csv(file1,check.names=F,row.names=1)
  data2 <- read.csv(file2,check.names=F,row.names=1)
  shared.parents <- setdiff(intersect(colnames(data1),colnames(data2)),
                            c("chrom","cM","bp"))
  ix <- which(data1$chrom==chrom)
  markers <- intersect(rownames(data1[ix,]), rownames(data2))
  data1 <- data1[markers, shared.parents]
  data2 <- data2[markers, shared.parents]
  np <- length(shared.parents)
  z <- strsplit(data1[,1],split="|",fixed=T)
  ploidy <- length(z[[1]])
  
  f <- function(x) {
    z <- strsplit(x,split="|",fixed=T)
    ploidy <- length(z[[1]])
    m <- length(z)
    y <- matrix(as.integer(unlist(z)),ncol=ploidy,nrow=m,byrow = T)
  }
  out <- NULL
  for (i in 1:np) {
    x1 <- f(data1[,i])
    colnames(x1) <- 1:ploidy
    x2 <- f(data2[,i])
    colnames(x2) <- ploidy + 1:ploidy
    x <- cbind(x1,x2)
    clustans <- hclust(dist(t(x)))
    cu <- cutree(clustans,k = ploidy)
    z <- split(as.integer(names(cu)),cu)
    z <- lapply(z,sort)
    x <- data.frame(parent=shared.parents[i], file1=sapply(z,"[[",1),
               file2=sapply(z,"[[",2))
    x$file2 <- x$file2 - ploidy
    out <- rbind(out,x)
    plot(as.dendrogram(clustans),horiz = TRUE, main=shared.parents[i],
        xlab="Genetic Distance")
  }
  rownames(out) <- NULL
  return(out)
}