#' Fit a single QTL model
#' 
#' Fit a single QTL model
#' 
#' The number of burn-in and total iterations in \code{params} can be estimated using \code{\link{set_params}}. Parameter \code{dominance} controls the genetic model for the QTL: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. The optional argument \code{cofactor} should be a list with three components: marker = name of the marker; dominance = 1, 2, 3, or 4; epistasis = TRUE/FALSE. When \code{polygenic = TRUE}, the model includes a random effect with covariance equal to the additive relationship computed by \code{\link{IBDmat}}, leaving out chromosome(s) with the QTL and cofactor (if present). Parameter \code{CI.prob} sets the probability (e.g., 0.90, 0.95) for the Bayesian credible interval for the estimated effects (to disable plotting of the CI, use \code{CI.prob=NULL}). 
#' 
#' The LOD and deltaDIC values returned by the function are relative to a model without \code{marker} but including the cofactor and polygenic effect when present. If \code{polygenic = FALSE}, the null model includes a GCA effect. r2 is the squared correlation between the fitted and observed values. The returned list \code{effects} contains the additive (and when included) digenic dominance effects. The proportion of variance for each effect is returned in \code{var}. The returned object \code{plots$dom} shows the digenic dominance effects above the diagonal, and below the diagonal is the sum of the additive and digenic dominance effects. 
#' 
#' @param data variable of class \code{\link{diallel_geno_pheno}}
#' @param trait name of trait
#' @param marker name of marker to fit as QTL
#' @param params list containing the number of burn-in (burnIn) and total iterations (nIter)
#' @param dominance dominance degree
#' @param cofactor optional, see Details for format.
#' @param polygenic TRUE/FALSE whether to include a polygenic effect
#' @param CI.prob probability for Bayesian credible interval
#' 
#' @return List containing
#' \describe{
#' \item{r2}{sauared correlation betwen fitted and observed values}
#' \item{deltaDIC}{Deviance Information Criterion relative to null model}
#' \item{resid}{Residuals}
#' \item{var}{Matrix with proportion of variance for the effects}
#' \item{effects}{List of matrices containing the additive and higher order effects}
#' \item{plots}{List of ggplot objects for the effects}
#' }
#' @examples
#' \dontrun{
#' ## additive effects
#' params1 <- set_params( diallel_example, trait = "tuber_shape" ,q=0.05,r=0.05)
#' 
#' fit1 <- fitQTL( data = diallel_example, 
#'                  trait = "tuber_shape", 
#'                  params = params1, 
#'                  marker = "solcap_snp_c2_25522",
#'                  CI.prob = 0.9)
#'                  
#' ## additive + dominance effects
#' params2 <- set_params( diallel_example, trait = "tuber_shape", dominance=2,q=0.05,r=0.05)
#' 
#' fit2 <- fitQTL( data = diallel_example, 
#'                  trait = "tuber_shape", 
#'                  params = params2, 
#'                  marker = "solcap_snp_c2_25522",
#'                  dominance = 2,
#'                  CI.prob=0.9)
#'                  
#' }
#'                  
#' @export
#' 
#' @importFrom stats var sd
#' @import ggplot2
#' @importFrom BGLR readBinMat

fitQTL <- function(data,trait,marker,params,dominance=1,cofactor=NULL,CI.prob=0.9,polygenic=TRUE) {
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  stopifnot(marker %in% data@map$marker)
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
  
  marker2 <- get_bin(marker,data@map)
  if (!is.null(cofactor)) {
    stopifnot(is.list(cofactor))
    stopifnot(cofactor$marker %in% data@map$marker)
    cofactor$marker <- get_bin(cofactor$marker,data@map)
    cofactor$X <- data@geno[[cofactor$marker]][1:cofactor$dominance]
    if (cofactor$epistasis) {
      cofactor$Xaa <- data@Z%*%faa(j=cofactor$marker,k=marker2,data=data)
    } else {
      cofactor$Xaa <- NULL
    }
    poly.marker <- c(marker,cofactor$marker)
  } else {
    poly.marker <- marker
  }
  
  if (polygenic) {
    chroms <- setdiff(unique(data@map$chrom),data@map$chrom[match(poly.marker,data@map$marker)])
    if (length(chroms)==0) {
      stop("There are no chromosomes remaining for the polygenic effect.")
    }
    G1 <- IBDmat(data=data,dominance=1,chrom=chroms)
    G1 <- data@Z %*% tcrossprod(G1,data@Z)
  } else {
    G1 <- NULL
  }
  
  y <- data@pheno[,trait]
  if (inherits(y,"factor")) {
    response <- "ordinal"
  } else {
    response <- "gaussian"
  }
  params <- list(response=response,nIter=params$nIter,burnIn=params$burnIn)

  #no marker model
  ans0 <- qtl1(y=y,X=data@X,Z=data@Z,params=params,X.GCA=data@X.GCA,cofactor=cofactor)

  #with marker
  ans1 <- qtl1(y=y,X=data@X,Z=data@Z,geno=data@geno[[marker2]][1:dominance],
               params=params,G1=G1,cofactor=cofactor)
  
  effect.lower <- effect.upper <- effect.mean <- vector("list",length=dominance)
  variances <- matrix(0,nrow=params$nIter-params$burnIn,ncol=dominance)
  colnames(variances) <- c("additive","digenic","trigenic","quadrigenic")[1:dominance]

  for (j in 1:dominance) {
    ans <- readBinMat(sub(pattern="X",replacement=j,x="tmp/ETA_aX_b.bin"))
    ans <- t(apply(ans,1,function(z){z-mean(z)})) #enforce constraint that sum of effects is zero
    effect.mean[[j]] <- apply(ans,2,mean)
    if (!is.null(CI.prob)) {
      tmp <- apply(ans,2,quantile,p=c(0.5-CI.prob/2,0.5+CI.prob/2))
      effect.lower[[j]] <- tmp[1,]
      effect.upper[[j]] <- tmp[2,]
    }
    variances[,j] <- apply(tcrossprod(data@Z %*% data@geno[[marker2]][[j]],ans),2,var)
  }
  
  if (polygenic) {
    var.poly <- scan("tmp/ETA_polyg_varU.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]*mean(diag(G1))
    variances <- cbind(variances,polygenic=var.poly)
  } 
  
  if (!is.null(cofactor)) {
    var.cof <- matrix(0,nrow=params$nIter-params$burnIn,ncol=cofactor$dominance)
    colnames(var.cof) <- paste("cofactor",c("additive","digenic","trigenic","quadrigenic")[1:cofactor$dominance],sep=".")
    for (j in 1:cofactor$dominance) {
      ans <- readBinMat(sub(pattern="X",replacement=j,x="tmp/ETA_cofX_b.bin"))
      var.cof[,j] <- apply(tcrossprod(data@Z %*% cofactor$X[[j]],ans),2,var)
    }
    variances <- cbind(variances,var.cof)
    if (cofactor$epistasis) {
      ans <- readBinMat("tmp/ETA_cofAA_b.bin")
      variances <- cbind(variances,cofactor.epistasis=apply(tcrossprod(cofactor$Xaa,ans),2,var))  
    }
  }
  
  varE <- scan("tmp/varE.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]
  h2 <- variances/(apply(variances,1,sum) + varE)  #proportion of variance for each term
  if (!is.null(CI.prob)) {
    tmp <- apply(h2,2,quantile,p=c(0.5-CI.prob/2,0.5+CI.prob/2))
    return.var <- cbind(Mean=round(apply(h2,2,mean),2),CI.lower=round(tmp[1,],2),CI.upper=round(tmp[2,],2))
  } else {
    return.var <- cbind(Mean=round(apply(h2,2,mean),2))
  }

  effects <- vector("list",length=ifelse(dominance==1,1,2))
  names(effects) <- c("add","dom")[1:length(effects)]

  if (!is.null(CI.prob)) {
    effects$add <- data.frame(Haplotype=attr(data@geno,"haplotypes"),Mean=effect.mean[[1]],CI.lower=effect.lower[[1]],CI.upper=effect.upper[[1]],stringsAsFactors = F)
  } else {
    effects$add <- data.frame(Haplotype=attr(data@geno,"haplotypes"),Mean=effect.mean[[1]],stringsAsFactors = F)
  }
  effects$add <- effects$add[order(effects$add$Haplotype),]

  if (dominance > 1) {
    diplotypes <- strsplit(attr(data@geno,"diplotypes"),split="+",fixed=T)
    if (!is.null(CI.prob)) {
      effects$dom <- data.frame(Haplotype1=sapply(diplotypes,function(x){x[1]}),Haplotype2=sapply(diplotypes,function(x){x[2]}),Mean=effect.mean[[2]],CI.lower=effect.lower[[2]],CI.upper=effect.upper[[2]],stringsAsFactors = F)
    } else {
      effects$dom <- data.frame(Haplotype1=sapply(diplotypes,function(x){x[1]}),Haplotype2=sapply(diplotypes,function(x){x[2]}),Mean=effect.mean[[2]],stringsAsFactors = F)
    }
    effects$dom <- effects$dom[order(effects$dom$Haplotype1,effects$dom$Haplotype2),]
  }

  #Plotting
  #Additive effects
  tmp <- strsplit(split=".",x=effects$add$Haplotype,fixed=T)
  founders <- sapply(tmp,function(x){x[1]})
  plotme <- data.frame(effects$add,founders=factor(founders))
  plotme$Haplotype <- rep(1:data@ploidy,length(levels(plotme$founders)))
  plotA <- ggplot(data=plotme,aes(x=Haplotype,y=Mean)) + 
    labs(title = paste("Trait:", trait),
         subtitle = paste("Marker:", marker)) + 
    theme_bw() + 
    scale_x_continuous(name="Haplotype",labels=1:data@ploidy,breaks=1:data@ploidy) +
    geom_bar(stat="identity",position="dodge",colour="black",fill="grey50") + 
    ylab("Additive Effect") + 
    theme(text = element_text(size=13)) + 
    facet_grid(.~founders,scales = "free_x")
  if (!is.null(CI.prob)) {
    plotA <- plotA + geom_errorbar(aes(ymax=CI.upper, ymin=CI.lower, width = 0.2))
  }
  
  if (dominance > 1) {
    #Digenic dominance effects
    plotme <- data.frame(x=c(effects$dom$Haplotype1,effects$dom$Haplotype2),y=c(effects$dom$Haplotype2,effects$dom$Haplotype1),z=c(effects$dom$Mean,effects$dom$Mean+effects$add$Mean[match(effects$dom$Haplotype1,effects$add$Haplotype)]+effects$add$Mean[match(effects$dom$Haplotype2,effects$add$Haplotype)])) #construct symmetric matrix
    plotme <- plotme[which(plotme[,1]!=plotme[,2]),] #remove diagonals
    plotD <- ggplot(data=plotme,aes(x=x,y=y,fill=z)) +
      geom_tile() +
      scale_fill_gradient2(name="") + 
      labs(title = paste("Trait:", trait),
         subtitle = paste("Marker:", marker)) +
      theme_bw() +
      xlab("") +
      ylab("") +
      theme(text = element_text(size=13),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
      coord_fixed(ratio=1)
      
    return(list(r2=round(ans1$r2,2),deltaDIC=ans1$DIC-ans0$DIC,resid=ans1$resid,var=return.var,effects=effects,plots=list(add=plotA,dom=plotD)))
  } else {
    return(list(r2=round(ans1$r2,2),deltaDIC=ans1$DIC-ans0$DIC,resid=ans1$resid,var=return.var,effects=effects,plots=list(add=plotA)))
  }
}

