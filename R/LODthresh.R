#' LOD thresholds for scan1
#' 
#' LOD thresholds for scan1
#' 
#' LOD thresholds to control the genome-wide false positive rate at 0.05 were determined via simulation for up to 20 parents and genome sizes up to 12 Morgans. A monotone increasing concave curve was fit to these results using R package \code{scam} and is used for prediction. (The LOD threshold does not depend on population size.)
#' 
#' @param genome.size Genome size in Morgans (not centiMorgans)
#' @param num.parents Number of parents
#' @param ploidy 2 or 4
#' @param dominance 1 (additive) or 2 (digenic dominance)
#' 
#' @return LOD threshold
#' @export
#' @importFrom scam predict.scam
#' 
LODthresh <- function(genome.size,num.parents,ploidy,dominance=1) {
  if (ploidy==2) {
    if (dominance==1) {
      ans <- predict.scam(diaQTL:::LOD2x.add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
    } else {
      ans <- predict.scam(diaQTL:::LOD2x.digenic,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
    }
  } 
  if (ploidy==4) {
    if (dominance==1) {
      ans <- predict.scam(diaQTL:::LOD4x.add,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
    } else {
      ans <- predict.scam(diaQTL:::LOD4x.digenic,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
    }
  }
  return(round(ans,2))
}