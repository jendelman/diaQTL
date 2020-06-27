#' Get parental allele dosage estimates
#' 
#' Get parental allele dosage estimates
#' 
#' Function can be used to get parental allele dosage estimates at a single marker for all individuals (in which case `id` should be `NULL`) or for a single individual for all markers (in which case `marker` should be `NULL`)
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param marker Name of marker
#' @param id Name of individual
#' 
#' @return Matrix of (id or markers) x parental alleles 
#' 
#' @examples
#' \dontrun{
#'   dosage_example = dosage(data = diallel_example, 
#'                           marker = "solcap_snp_c2_25522")
#'   dosage_example = dosage(data = diallel_example, 
#'                           id = "W15263-8R")
#' }
#' 
#' @export

dosage <- function(data,marker=NULL,id=NULL) {
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(is.null(marker)|is.null(id))
  stopifnot(!is.null(marker)|!is.null(id))
  if (!is.null(marker)) {
    ans <- data@geno$A[[marker]]
    rownames(ans) <- attr(data@geno$A,"id")
    colnames(ans) <- attr(data@geno$A,"alleles")
    return(ans)
  }
  if (!is.null(id)) {
    k <- match(id,attr(data@geno$A,"id"))
    ans <- t(sapply(data@geno$A,function(geno){geno[k,]}))
    rownames(ans) <- names(data@geno$A)
    colnames(ans) <- attr(data@geno$A,"alleles")
    return(ans)
  }
}
