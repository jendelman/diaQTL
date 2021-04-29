#' Generate diaQTL input files from OneMap
#' 
#' Generate diaQTL input files from 'onemap_progeny_haplotypes' object class of OneMap R package (version >2.2.0)
#' 
#' @param onemap.object onemap_progeny_haplotypes object class
#' @param outstem prefix for the pedigree and genotype files for diaQTL
#' @param conv.bp value to multiply cM of the map position to have in base pairs
#' 
#' @return NULL
#' 
#' #' @examples
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
                           conv.bp=10^6,
                           outstem=""){
  
  if(!("onemap_progeny_haplotypes" %in% class(data)))
    stop("data should be from 'onemap_progeny_haplotypes' class. Use 'progeny_haplotypes' function to create such object (OneMap > v2.2.0)")
  
  data$haplotypes = paste0(data$parents,data$homologs)
  data = dcast(data, ind + marker + grp + pos ~ haplotypes, value.var="prob")
  data$G1 = round(data$P1H1*data$P2H1,digits)
  data$G2 = round(data$P1H1*data$P2H2,digits)
  data$G3 = round(data$P1H2*data$P2H1,digits)
  data$G4 = round(data$P1H2*data$P2H2,digits)
  
  ## removing the 0s at the end
  check = data$G1+data$G2+data$G3+data$G4
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
                    bp=round(data$pos,digits)*conv.bp,
                    prob=data$prob)
  data = dcast(data, marker + chrom + cM + bp ~ ind, value.var="prob")
  data = data[order(data$chrom,data$cM),]
  write.csv(data,file=paste0(outstem,"diaQTL_geno.csv"),row.names=F)
  
  ped = data.frame(id=colnames(data)[-c(1:4)],
                   parent1="P1",
                   parent2="P2")
  write.csv(ped,file=paste0(outstem,"diaQTL_ped.csv"),row.names=F)
}