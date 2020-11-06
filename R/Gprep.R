#' Compute realized relationship matrix for fitQTL
#' 
#' Compute realized relationship matrix for fitQTL
#' 
#' Computes the additive realized relationship matrix with \code{\link{IBDmat}} so that \code{\link{fitQTL}} does not need to and can run more quickly. The chromosome containing \code{marker} is excluded from the relationship calculation.
#' 
#' @param data variable inheriting from class \code{\link{diallel_geno_pheno}}
#' @param marker marker that will be used for \code{\link{fitQTL}} 
#' 
#' @return variable inheriting from class \code{\link{diallel_geno_pheno_G}}
#' 
#' @examples
#' \dontrun{
#'   diallel_example <- Gprep(data = diallel_example, marker = "solcap_snp_c2_25522")
#' }
#' 
#' @export

Gprep <- function(data,marker) {
  stopifnot(marker %in% data@map$marker)
  stopifnot(inherits(data,"diallel_geno_pheno"))
  G1 <- IBDmat(data=data,dominance=1,chrom=setdiff(unique(data@map$chrom),data@map$chrom[match(marker,data@map$marker)]))
  return(new(Class="diallel_geno_pheno_G",data,G1=G1))
}