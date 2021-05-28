#' delta DIC thresholds for scan1
#' 
#' delta DIC thresholds for scan1
#' 
#' Thresholds to control the genome-wide false positive rate at `alpha` were determined for half-diallel mating designs with up to 10 parents. 
#' 
#' @param genome.size Genome size in Morgans (not centiMorgans)
#' @param num.parents Number of parents (2 to 10)
#' @param ploidy 2 or 4
#' @param alpha significance level (0.01, 0.05, 0.10, or 0.20)
#' @param dominance 1 (additive) or 2 (digenic dominance)
#' 
#' @return -deltaDIC threshold
#' @examples
#' \dontrun{
#'   DIC_thresh(genome.size=10, 
#'              num.parents=4,
#'              ploidy=4,
#'              dominance=1,
#'              alpha=0.05)
#'   } 
#'   
#' @export
#' @importFrom scam predict.scam
#' 

DIC_thresh <- function(genome.size,num.parents,ploidy,alpha=0.05,dominance=1){

  ## Checking for input errors
  if(!is.element(num.parents,2:10))
    stop(deparse("number of parents should be between 2 and 10"))
  
  if(!c(ploidy %in% c(2,4)))
    stop(deparse("ploidy should be 2 or 4"))
  
  if(!c(alpha %in% c(0.01,0.05,0.10,0.20)))
    stop(deparse("alpha should be 0.01, 0.05, 0.10, or 0.20"))
  
  if(!c(dominance %in% c(1,2)))
    stop(deparse("dominance should be 1 (additive) or 2 (additive+digenic dominance)"))
  
  if(alpha==0.01){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x01add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x01dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x01add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x01dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
  }
  if(alpha==0.05){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x05add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x05dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x05add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x05dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
  }
  
  if(alpha==0.10){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x10add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x10dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x10add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x10dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
  }
  
  if(alpha==0.20){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x20add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x20dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x20add,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x20dig,newdata=data.frame(Genome.Size=log10(genome.size),Parents=num.parents))
      }
    }
  }
  
  return(-round(ans,2))
}
