#' Generate diaQTL input files from MAPpoly
#' 
#' Generate diaQTL input files from a list of 'mappoly.genoprob' object class of 
#' MAPpoly R package (version >0.2.3). It only works for autotetraploid data.
#' 
#' @param MAPpoly.object list of mappoly.genoprob object class (one for each linkage group)
#' @param outstem prefix for the pedigree and genotype files for diaQTL
#' @param conv.bp value to multiply cM of the map position to have in base pairs
#' 
#' @return NULL
#' 
#' #' @examples
#' \dontrun{
#'     # see MAPpoly tutorial for details on its functions
#'     MAP <- list(lg1, lg2, lg3) #your map with your linkage groups
#'     genoprob <- vector("list", 3) #if 3 linkage groups
#'     for(i in 1:length(genoprob))
#'         genoprob[[i]] <- mappoly::calc_genoprob_error(input.map = MAP[[i]], error = 0.05)
#'     convert_mappoly(genoprob)
#' }
#' @export
#' @importFrom utils write.csv
#' @importFrom reshape2 dcast
#' @importFrom plyr adply
#' 


convert_mappoly <- function(data,
                            digits=4,
                            conv.bp=10^6,
                            outstem=""){
  
  if(unique(unlist(lapply(data,class)))!="mappoly.genoprob")
    stop("data should be a list of 'mappoly.genoprob' class. Use 'calc_genoprob_error' function to create such object (MAPpoly > v0.2.3.1)")
  
  outputAll = NULL  
  for(i in 1:length(data)){
    output = NULL
    output = data.frame(marker=names(data[[i]]$map),
                        chrom=i,
                        cM=round(data[[i]]$map,digits),
                        bp=round(data[[i]]$map,digits)*conv.bp)
    data[[i]]$probs=round(data[[i]]$probs,4)
    tmp = adply(data[[i]]$probs,c(2,3))
    if(i==1){ #get the genotypic code from diaQTL 
      genolabels = names(tmp)[-c(1:2)]          
      for(j in 1:8)
        genolabels = gsub(letters[j],j,genolabels)
      genolabels = gsub(":","",genolabels)
      genolabels = match(genolabels,gsub("-","",F1codes$State))
      genolabels = paste0(genolabels,collapse="|")
    }
    
    for(j in 3:37)
      tmp$prob = paste0(tmp$prob,format(tmp[,j],scientific=FALSE),"|")
    tmp$prob = paste0("=>",tmp$prob,format(tmp[,38],scientific=FALSE))
    tmp$prob = paste0(genolabels,tmp$prob) #paste string
    
    ## Trimming zeros
    tmp$prob = gsub(paste0(c("0.",paste0(rep(0,digits),collapse="")),collapse=""),
                    "0",
                    tmp$prob)
    tmp = tmp[,-c(3:38)]
    names(tmp) = c("marker","ind","prob")
    tmp = dcast(tmp, marker ~ ind, value.var="prob")
    output = cbind(output,tmp[,-1])
    
    outputAll = rbind(outputAll,output)
  }
  
  write.csv(outputAll,file=paste0(outstem,"diaQTL_geno.csv"),row.names=F)
  
  ped = data.frame(id=colnames(outputAll)[-c(1:4)],
                   parent1="P1",
                   parent2="P2")
  write.csv(ped,file=paste0(outstem,"diaQTL_ped.csv"),row.names=F)
}
