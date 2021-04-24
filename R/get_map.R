#' Get map summary from diallel_geno object
#' 
#' Get map summary from diallel_geno object
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param summary logical, if TRUE (default) returns total sizes per chromosome, if FALSE returns the map
#' 
#' @return data frame with map summary or the map
#' 
#' @examples
#' \dontrun{
#'   get_map(diallel_example)
#' }
#' 
#' @export
#' 

get_map <- function(data,summary=TRUE){
  if(!(class(data) %in% c("diallel_geno","diallel_geno_pheno")))
    stop("data should be diallel_geno or diallel_geno_pheno")
  
  if(summary){
    tmp = tapply(data@map$cM,data@map$chrom,max)
    total = sum(tmp)
    return(data.frame(chrom=c(names(tmp),"total:"),
                      'Morgans'=round((c(as.numeric(tmp),total)/100),2)))
  }else{
    return(data@map)
  }
}