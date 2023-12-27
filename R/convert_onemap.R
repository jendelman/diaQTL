#' Generate diaQTL input files from OneMap
#' 
#' Generate diaQTL input files from 'onemap_progeny_haplotypes' object class of OneMap R package (version >2.2.0)
#' 
#' @param data onemap_progeny_haplotypes object class
#' @param digits how many rounding digits for the probabilities output (default=4)
#' @param outstem prefix for the pedigree and genotype files for diaQTL
#' 
#' @return NULL
#' 
#' @examples
#' \dontrun{
#'     map <- list(LG1_final, LG2_final)
#'     progeny_haplot <- onemap::progeny_haplotypes(map,
#'                                                  most_likely = FALSE,
#'                                                  ind = "all")
#'     convert_onemap(progeny_haplot)
#' }
#' @export
#' @importFrom utils write.csv
#' @importFrom reshape2 dcast

convert_onemap <- function(data,
                           digits=4,
                           outstem=""){
  
  if(!("onemap_progeny_haplotypes" %in% class(data)))
    stop("data should be from 'onemap_progeny_haplotypes' class. Use 'progeny_haplotypes' function to create such object (OneMap > v2.2.0)")
  
  data$haplotypes = paste0(data$parents,data$parents.homologs)
  data = dcast(data, ind + marker + grp + pos ~ haplotypes, value.var="prob")
  data$G1 = round(data$P1H1*data$P2H1,digits)
  data$G2 = round(data$P1H1*data$P2H2,digits)
  data$G3 = round(data$P1H2*data$P2H1,digits)
  data$G4 = round(data$P1H2*data$P2H2,digits)
  
  ## removing the 0s at the end
  check = data$G1+data$G2+data$G3+data$G4
  if(length(which(check==0)) > 0)
    data = data[-which(check==0),]
  
  data$prob = paste0("1|2|3|4=>",
                     format(data$G1,scientific=FALSE),"|",
                     format(data$G2,scientific=FALSE),"|",
                     format(data$G3,scientifi=FALSE),"|",
                     format(data$G4,scientific=FALSE))
  data = data.frame(ind=data$ind,
                    marker=data$marker,
                    chrom=data$grp,
                    cM=round(data$pos,digits),
                    prob=data$prob)
  data = dcast(data, marker + chrom + cM ~ ind, value.var="prob")
  data = data[order(data$chrom,data$cM),]
  write.csv(data,file=paste0(outstem,"diaQTL_geno.csv"),row.names=F)
  
  ped = data.frame(id=colnames(data)[-c(1:3)],
                   parent1="P1",
                   parent2="P2")
  write.csv(ped,file=paste0(outstem,"diaQTL_ped.csv"),row.names=F)
}