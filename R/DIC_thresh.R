#' delta DIC thresholds for scan1
#' 
#' delta DIC thresholds for scan1
#' 
#' delta DIC thresholds to control the genome-wide false positive rate at level `alpha` were determined via simulation for up to 20 parents and genome sizes up to 12 Morgans. A monotone decreasing concave curve was fit to these results using R package \code{scam} and is used for prediction. 
#' 
#' @param genome.size Genome size in Morgans (not centiMorgans)
#' @param num.parents Number of parents
#' @param ploidy 2 or 4
#' @param alpha false positive rate: 0.01, 0.05, 0.10, or 0.20
#' @param dominance 1 (additive) or 2 (digenic dominance)
#' 
#' @return deltaDIC threshold
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
  if(!is.numeric(genome.size) | genome.size>12.5)
    stop(deparse("genome size should be numeric (in Morgans) and lower than 12"))
  
  if(!is.numeric(num.parents) | (round(num.parents)!=num.parents))
    stop(deparse("number of parents should be integer and numeric"))
  
  if(!c(ploidy %in% c(2,4)))
    stop(deparse("ploidy should be 2 or 4"))
  
  if(!c(alpha %in% c(0.01,0.05,0.10,0.20)))
    stop(deparse("alpha should be 0.01, 0.05, 0.10, or 0.20"))
  
  if(!c(dominance %in% c(1,2)))
    stop(deparse("dominance should be 1 (additive) or 2 (additive+digenic dominance)"))
  
  if(alpha==0.01){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x01add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x01dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x01add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x01dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
  }
  if(alpha==0.05){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x05add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x05dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x05add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x05dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
  }
  
  if(alpha==0.10){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x10add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x10dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x10add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x10dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
  }
  
  if(alpha==0.20){
    if (ploidy==2) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC2x20add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC2x20dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
    if (ploidy==4) {
      if (dominance==1) {
        ans <- predict.scam(deltaDIC4x20add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      } else {
        ans <- predict.scam(deltaDIC4x20dig,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
      }
    }
  }
  
  return(round(-ans,2))
}