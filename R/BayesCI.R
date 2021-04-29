#' Bayesian Credible Interval for QTL position
#' 
#' Bayesian Credible Interval for QTL position
#' 
#' Parameter \code{CI.prob} sets the probability for the Bayesian credible interval (e.g., 0.90, 0.95) 
#' using the profile likelihood (posterior mean).
#' 
#' @param scan1_data data frame output from scan1
#' @param data variable of class \code{\link{diallel_geno_pheno}}
#' @param chrom chromosome
#' @param CI.prob probability for the credible interval
#' 
#' @return subset of scan1_data with markers in the CI
#'
#' @examples
#' \dontrun{
#'   BayesCI(scan1_example, diallel_example, chrom="10", CI.prob=0.9)
#'   }
#' @export

BayesCI <- function(scan1_data,data,chrom,CI.prob=0.9) {
  stopifnot(chrom %in% scan1_data$chrom)
  
  scan1a <- scan1_data[match(names(data@geno),scan1_data$marker),]
  ix <- which(scan1a$chrom == chrom)
  x <- scan1a$cM[ix]

  tmp <- data.frame(pos=ix,left=c(0,diff(x)/2),right=c(diff(x)/2,0),prob=exp(scan1a$LL[ix]))

  tmp$y <- tmp$prob*(tmp$right+tmp$left)
  tmp$area <- tmp$y/sum(tmp$y)
  n <- length(x)
  tmp$cdf <- apply(array(1:n),1,function(k){sum(tmp$area[1:k])})
  
  lower <- max(which(tmp$cdf <= 0.5-CI.prob/2))
  lower <- max(1,lower-1)
  upper <- min(which(tmp$cdf >= 0.5+CI.prob/2))
  upper <- min(upper+1,n)
  limits <- scan1a$cM[ix[c(lower,upper)]]
  out <- scan1_data[scan1_data$chrom == chrom & scan1_data$cM >= limits[1] & scan1_data$cM <= limits[2],]
  return(out)
}
