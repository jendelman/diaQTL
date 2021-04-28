DIC_thresh <- function(genome.size,num.parents,ploidy,alpha=0.05,dominance=1){

  ## Checking for input errors  
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
  
  return(round(ans,2))
}
