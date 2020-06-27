#' LOD thresholds for scan1
#' 
#' LOD thresholds for scan1
#' 
#' LOD thresholds to control the genome-wide false positive rate at 0.05 were determined via simulation for up to 20 parents and genome sizes up to 12 Morgans. A monotone increasing concave curve was fit to these results using R package \code{scam} and is used for prediction. (The LOD threshold does not depend on population size.)
#' 
#' @param genome.size Genome size in Morgans (not centiMorgans)
#' @param num.parents Number of parents
#' @param ploidy 2 or 4
#' 
#' @return LOD threshold
#' @export
#' @importFrom scam predict.scam
#' 
LODthresh <- function(genome.size,num.parents,ploidy) {
  if (ploidy==2) {
    ans <- predict(diaQTL:::LOD2x,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
  } 
  if (ploidy==4) {
    ans <- predict(diaQTL:::LOD4x,newdata=data.frame(Genome.Size=genome.size,Num.Parents=num.parents))
  }
  return(round(ans,1))
}