#' Generate diaQTL input files from MAPpoly
#' 
#' Generate diaQTL input files from a list of 'mappoly.genoprob' object class of 
#' MAPpoly R package (version >0.2.3). It only works for autotetraploid data.
#' 
#' @param data list of mappoly.genoprob object class (one for each linkage group)
#' @param ploidy Either 2 or 4 (default=4)
#' @param digits how many rounding digits for the probabilities output (default=4)
#' @param outstem prefix for the pedigree and genotype files for diaQTL
#' 
#' @return NULL
#' 
#' @examples
#' \dontrun{
#'     # see MAPpoly tutorial for details on its functions
#'     MAP <- list(lg1, lg2, lg3) #your map with your linkage groups
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

convert_mappoly <- function(data,
                            digits=4,
                            ploidy=4,
                            outstem=""){
  
  if(unique(unlist(lapply(data,class)))!="mappoly.genoprob")
    stop("data should be a list of 'mappoly.genoprob' class. Use 'calc_genoprob_error' function to create such object (MAPpoly > v0.2.3.1)")
  
  outputAll = NULL  
  if(ploidy==4){
    for(i in 1:length(data)){
      output = NULL
      output = data.frame(marker=names(data[[i]]$map),
                          chrom=i,
                          cM=round(data[[i]]$map,digits))
      data[[i]]$probs=round(data[[i]]$probs,digits)
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
  }else{
    F1codesDiploid = c("1-3","1-4","2-3","2-4")
    for(i in 1:length(data)){
      output = NULL
      output = data.frame(marker=names(data[[i]]$map),
                          chrom=i,
                          cM=round(data[[i]]$map,digits))
      data[[i]]$probs=round(data[[i]]$probs,digits)
      tmp = adply(data[[i]]$probs,c(2,3))
      if(i==1){ #get the genotypic code from diaQTL 
        genolabels = names(tmp)[-c(1:2)]          
        for(j in 1:4)
          genolabels = gsub(letters[j],j,genolabels)
        genolabels = gsub(":","",genolabels)
        genolabels = match(genolabels,gsub("-","",F1codesDiploid))
        genolabels = paste0(genolabels,collapse="|")
      }
      
      for(j in 3:5)
        tmp$prob = paste0(tmp$prob,format(tmp[,j],scientific=FALSE),"|")
      tmp$prob = paste0("=>",tmp$prob,format(tmp[,6],scientific=FALSE))
      tmp$prob = paste0(genolabels,tmp$prob) #paste string
      
      ## Trimming zeros
      tmp$prob = gsub(paste0(c("0.",paste0(rep(0,digits),collapse="")),collapse=""),
                      "0",
                      tmp$prob)
      tmp = tmp[,-c(3:6)]
      names(tmp) = c("marker","ind","prob")
      tmp = dcast(tmp, marker ~ ind, value.var="prob")
      output = cbind(output,tmp[,-1])
      
      outputAll = rbind(outputAll,output)
    }
  }
    
  write.csv(outputAll,file=paste0(outstem,"diaQTL_geno.csv"),row.names=F)
  
  ped = data.frame(id=colnames(outputAll)[-c(1:3)],
                   parent1="P1",
                   parent2="P2")
  write.csv(ped,file=paste0(outstem,"diaQTL_ped.csv"),row.names=F)
}
