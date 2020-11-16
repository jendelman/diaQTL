#' Fit a single QTL model
#' 
#' Fit a single QTL model
#' 
#' For quantitative traits, R2 is the percent of variation explained by the regression (MSS/TSS). For binary traits, R2 is the squared phi correlation (as a percentage). LOD score is the difference between the log10-likelihood for the QTL model vs. no QTL model; higher values are better. deltaDIC is the difference between the Deviance Information Criterion for the QTL model vs. no QTL model; lower values are better. Parameter \code{dominance} controls the genetic model: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. MCMC \code{params} can be estimated using \code{\link{set_params}}. Parameter \code{CI.prob} sets the probability (e.g., 0.90, 0.95) for the Bayesian credible interval for the estimated effects. The returned list \code{effects} contains the additive (and when included) digenic dominance effects. The proportion of variance for each effect is returned in \code{var}. The returned object \code{plots$dom} shows the digenic dominance effects above the diagonal, and below the diagonal is the sum of the additive and digenic dominance effects. The polygenic background effect has covariance equal to the additive relationship computed by \code{\link{IBDmat}}, leaving out the chromosome with the QTL. For faster execution with the polygenic model, use \code{\link{Gprep}} first.
#' 
#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param trait Name of trait
#' @param marker Name of marker to fit as QTL
#' @param params List containing the number of burn-in (burnIn) and total iterations (nIter)
#' @param dominance Dominance degree 
#' @param CI.prob Probability for Bayesian credible interval
#' @param polygenic TRUE/FALSE whether to include polygenic background effect
#' @param cofactor Name of marker to fit as cofactor (optional)
#' 
#' @return List containing
#' \describe{
#' \item{R2}{Coefficient of determination}
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

fitQTL <- function(data,trait,marker,params,dominance=1,CI.prob=0.9,polygenic=TRUE,cofactor=NULL) {

  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  stopifnot(marker %in% data@map$marker)
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
  if (polygenic) {
    if (!inherits(data,"diallel_geno_pheno_G")) {
      data <- Gprep(data,marker)
    }
    G1 <- data@Z %*% tcrossprod(data@G1,data@Z)
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

  if (!is.null(cofactor)) {
    stopifnot(cofactor %in% data@map$marker)
    Xcof <- data@Z%*%data@geno[[get_bin(cofactor,data@map)]][[1]]
  } else {
    Xcof <- NULL
  }
  
  #no marker model
  ans0 <- qtl1(y=y,X=data@X,Z=data@Z,params=params,Xcof=Xcof,X.GCA=data@X.GCA)
  
  #with marker
  marker2 <- get_bin(marker,data@map)
  ans1 <- qtl1(y=y,X=data@X,Z=data@Z,geno=data@geno[[marker2]][1:dominance],params=params,Xcof=Xcof,G1=G1)
  
  effect.lower <- effect.upper <- effect.mean <- vector("list",length=dominance)
  variances <- matrix(0,nrow=params$nIter-params$burnIn,ncol=dominance+as.integer(polygenic))
  colnames(variances) <- c("polygenic","additive","digenic","trigenic","quadrigenic")[1:(dominance+as.integer(polygenic))]
  
  if (polygenic) {
    variances[,1] <- scan("tmp/ETA_polyg_varU.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]*mean(diag(G1))
  } 
  
  for (j in 1:dominance) {
    ans <- readBinMat(sub(pattern="X",replacement=j,x="tmp/ETA_aX_b.bin"))
    effect.mean[[j]] <- apply(ans,2,mean)
    tmp <- apply(ans,2,quantile,p=c(0.5-CI.prob/2,0.5+CI.prob/2))
    effect.lower[[j]] <- tmp[1,]
    effect.upper[[j]] <- tmp[2,]
    variances[,as.integer(polygenic)+j] <- apply(tcrossprod(data@Z %*% data@geno[[marker2]][[j]],ans),2,var)
  }
  varE <- scan("tmp/varE.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]
  h2 <- variances/(apply(variances,1,sum) + varE)  #proportion of variance for each term
  tmp <- apply(h2,2,quantile,p=c(0.5-CI.prob/2,0.5+CI.prob/2))
  return.var <- cbind(Mean=round(apply(h2,2,mean),2),CI.lower=round(tmp[1,],2),CI.upper=round(tmp[2,],2))

  effects <- vector("list",length=ifelse(dominance==1,1,2))
  names(effects) <- c("add","dom")[1:length(effects)]

  effects$add <- data.frame(Haplotype=attr(data@geno,"haplotypes"),Mean=effect.mean[[1]],CI.lower=effect.lower[[1]],CI.upper=effect.upper[[1]],stringsAsFactors = F)
  effects$add <- effects$add[order(effects$add$Haplotype),]

  if (dominance > 1) {
    diplotypes <- strsplit(attr(data@geno,"diplotypes"),split="+",fixed=T)
    effects$dom <- data.frame(Haplotype1=sapply(diplotypes,function(x){x[1]}),Haplotype2=sapply(diplotypes,function(x){x[2]}),Mean=effect.mean[[2]],CI.lower=effect.lower[[2]],CI.upper=effect.upper[[2]],stringsAsFactors = F)
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
    geom_bar(stat="identity",position="dodge",colour="black") + 
    ylab("Additive Effect") + 
    theme(text = element_text(size=13)) + 
    facet_grid(.~founders,scales = "free_x") +
    geom_errorbar(aes(ymax=CI.upper, ymin=CI.lower, width = 0.2))
  
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
      
    return(list(R2=round(ans1$R2,2),deltaDIC=ans1$DIC-ans0$DIC,resid=ans1$resid,var=return.var,effects=effects,plots=list(add=plotA,dom=plotD)))
  } else {
    return(list(R2=round(ans1$R2,2),deltaDIC=ans1$DIC-ans0$DIC,resid=ans1$resid,var=return.var,effects=effects,plots=list(add=plotA)))
  }
}

