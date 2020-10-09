#' Dosage of parental diplotypes
#' 
#' Dosage of parental diplotypes
#' 
#' Function can be used to get parental diplotype dosage estimates at a single marker for all individuals (in which case `id` should be `NULL`) or for a single individual for all markers (in which case `marker` should be `NULL`)
#' 
#' @param data Variable inheriting from class \code{\link{diallel_geno}}
#' @param marker Name of marker
#' @param id Name of individual
#' 
#' @return Matrix of (id or markers) x parental diplotypes 
#' 
#' @examples
#' \dontrun{
#'   diplo_example = diplo_get(data = diallel_example, 
#'                           marker = "solcap_snp_c2_25522")
#'   diplo_example = diplo_get(data = diallel_example, 
#'                           id = "W15263-8R")
#' }
#' 
#' @export

diplo_get <- function(data,marker=NULL,id=NULL) {
  stopifnot(inherits(data,"diallel_geno"))
  stopifnot(is.null(marker)|is.null(id))
  stopifnot(!is.null(marker)|!is.null(id))
  if (!is.null(marker)) {
    k <- get_bin(marker,data@map)
    if (is.na(k)) {stop("Marker not present")}
    ans <- data@geno[[k]][[2]]
    rownames(ans) <- attr(data@geno,"id")
    colnames(ans) <- attr(data@geno,"diplotypes")
    return(as.matrix(ans))
  }
  if (!is.null(id)) {
    k <- match(id,attr(data@geno,"id"))
    if (is.na(k)) {stop("Individual not present")}
    ans <- t(sapply(data@geno,function(geno){geno[[1]][k,]}))
    colnames(ans) <- attr(data@geno,"diplotypes")
    ans <- ans[data@map$bin,] #expand to all markers
    rownames(ans) <- data@map$marker
    return(as.matrix(ans))
  }
}
