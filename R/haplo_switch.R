#' Find haplotype switches
#' 
#' Find haplotype switches
#' 
#' Designed to help with fine mapping of QTL. Function returns the location of the nearest haplotype switch on both sides (left, right) of a marker and haplotype of interest that represent the presumed location of a QTL allele.  
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param marker Name of focal marker 
#' @param haplotype Name of parental haplotype
#' @param position Either "cM" or "bp"
#' @param jump Change in dosage to designate haplotype switch
#' 
#' @return Data frame with the locations of the nearest haplotype switch on the left and right side of the focal marker for all individuals with a haplotype switch
#' 
#' @export

haplo_switch <- function(data,marker,haplotype,position,jump=0.8) {
  
  stopifnot(haplotype %in% attr(data@geno,"haplotypes"))
  stopifnot(marker %in% data@map$marker)
  stopifnot(position %in% colnames(data@map))
  
  id <- attr(data@geno,"id")
  n <- length(id)
  
  marker <- get_bin(marker,data@map)
  chrom <- data@map$chrom[match(marker,data@map$marker)]
  map <- data@map[data@map$chrom==chrom,]
  x0 <- map[match(marker,map$marker),position]
  result <- data.frame(id=id,dosage=numeric(n),left.marker=as.character(rep(NA,n)),left.pos=as.numeric(rep(NA,n)),right.marker=as.character(rep(NA,n)),right.pos=as.numeric(rep(NA,n)),stringsAsFactors = F)
  for (j in 1:n) {
    y <- haplo_get(data=data,id=id[j])
    y <- y[map$marker,haplotype]
    y2 <- abs(y-y[marker])
    result$dosage[j] <- y[marker]
    iq <- which(map[,position] > x0 & y2 > jump)
    if (length(iq)>0) {
      result$right.pos[j] <- min(map[iq,position])
      result$right.marker[j] <- map$marker[iq[which.min(map[iq,position])]]
    }
    iq <- which(map[,position] < x0 & y2 > jump)
    if (length(iq)>0) {
      result$left.pos[j] <- max(map[iq,position])
      result$left.marker[j] <- map$marker[iq[which.max(map[iq,position])]]
    }
  }
  return(result[which(!is.na(result$left.marker) | !is.na(result$right.marker)),])
}
