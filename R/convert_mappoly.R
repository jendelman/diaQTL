#' Generate diaQTL input files from MAPpoly
#' 
#' Generates diaQTL input files from the output of calc_genoprob or calc_genoprob_error in the MAPpoly package. The argument \code{data} is a list containing the results for each linkage group. Map distances are rounded to 0.01 cM, and genotype probabilties are rounded to three decimal places. 
#' 
#' @param data list of variables of class mappoly.genoprob (one for each linkage group)
#' @param ploidy Either 2 or 4 
#' @param outstem prefix for the pedigree and genotype files for diaQTL
#' 
#' @return NULL
#' 
#' @examples
#' \dontrun{
#'     # see MAPpoly tutorial for details on its functions
#'     MAP <- list(lg1, lg2, lg3) #list of linkage groups
#'     genoprob <- vector("list", 3) #if 3 linkage groups
#'     for(i in 1:length(genoprob))
#'         genoprob[[i]] <- mappoly::calc_genoprob_error(input.map = MAP[[i]], error = 0.05)
#'     convert_mappoly(genoprob, ploidy=4)
#' }
#' @export
#' @importFrom utils write.csv
#' @importFrom reshape2 dcast
#' @importFrom plyr adply
#' 

convert_mappoly <- function(data, ploidy, outstem=""){
  
  stopifnot(ploidy %in% c(2,4))
  #if(unique(unlist(lapply(data,class)))!="mappoly.genoprob")
  #  stop("data should be a list of 'mappoly.genoprob' class")
  
  genolabels <- dimnames(data[[1]]$probs)[[1]]
  for (j in 1:(2*ploidy))
    genolabels <- gsub(letters[j],j,genolabels)
  genolabels <- gsub(":","",genolabels)
  
  if (ploidy==4) {
    states <- gsub("-","",F1codes$State)
  } else {
    states <- c("13","14","23","24")
  }
  
  outputAll <- NULL 
  for (i in 1:length(data)) {
    tmp <- adply(round(data[[i]]$probs,3),c(2,3))
    colnames(tmp) <- c("marker","id",match(genolabels,states))
    prob <- apply(tmp[,-(1:2)],1,function(x){x[x>0]})
    tmp$prob <- sapply(prob,function(x){
      paste(paste0(names(x),collapse="|"),paste0(x,collapse="|"),sep="=>")
      })
    tmp <- dcast(tmp[,c("marker","id","prob")], marker ~ id, value.var="prob")
    output <- data.frame(marker=names(data[[i]]$map),
                        chrom=i,
                        cM=round(data[[i]]$map,2),tmp[,-1])
    
    outputAll = rbind(outputAll,output)
  }

  write.csv(outputAll,file=paste0(outstem,"diaQTL_geno.csv"),row.names=F)
  
  ped <- data.frame(id=colnames(outputAll)[-c(1:3)],
                   parent1="P1",
                   parent2="P2")
  write.csv(ped,file=paste0(outstem,"diaQTL_ped.csv"),row.names=F)
}
