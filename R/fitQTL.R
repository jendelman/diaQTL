#' Fit a single QTL model
#' 
#' Fit a single QTL model
#' 
#' Standard errors of the posterior mean estimates are calculated by dividing the SD of the Markov Chain by the square root of the effective number of iterations, which is calculated by function \code{effectiveSize} in R package \code{coda}. The error bars on the plot of additive effects correspond to +/- 1.96*SE (95 percent confidence interval). For binary traits, R2 = the squared phi correlation. Parameter \code{dominance} controls the genetic model: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. 
#' 
#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param params List containing the number of burn-in (burnIn) and total iterations (nIter)
#' @param dominance Dominance degree (1-4). See Details.
#' @param trait Name of trait
#' @param marker Name of marker to fit as QTL
#' @param cofactor Name of marker to fit as cofactor
#' 
#' @return List containing
#' \describe{
#' \item{R2}{Coefficient of determination}
#' \item{deltaDIC}{Deviance Information Criterion relative to null model}
#' \item{var}{Matrix with proportion of variance for additive (= heritability) and higher order effects}
#' \item{effects}{List of matrices containing the additive and higher order effects}
#' \item{plots}{List of ggplot objects for the effects}
#' }
#' @examples
#' \dontrun{
#' ## additive effects
#' params1 <- set_params( diallel_example, trait = "tuber_shape" )
#' 
#' fit1 <- fitQTL( data = diallel_example, 
#'                  trait = "tuber_shape", 
#'                  params = params1, 
#'                  marker = "solcap_snp_c2_25522")
#'                  
#' ## additive + dominance effects
#' params2 <- set_params( diallel_example, trait = "tuber_shape", dominance = TRUE )
#' 
#' fit2 <- fitQTL( data = diallel_example, 
#'                  trait = "tuber_shape", 
#'                  params = params2, 
#'                  marker = "solcap_snp_c2_25522",
#'                  dominance = 2)
#'                  
#' }
#'                  
#' @export
#' 
#' @importFrom stats var sd
#' @importFrom coda effectiveSize
#' @import ggplot2
#' @importFrom BGLR readBinMat

fitQTL <- function(data,params,dominance=1,trait,marker,cofactor=NULL) {
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  stopifnot(marker %in% data@map$marker)
  if (dominance > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
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

  k <- get_bin(marker,data@map)
  #with marker
  ans1 <- qtl1(y=y,X=data@X,Z=data@Z,geno=data@geno[[k]][1:dominance],params=params,Xcof=Xcof)
  
  effect.SE <- effect.mean <- vector("list",length=dominance)
  variances <- matrix(0,nrow=params$nIter-params$burnIn,ncol=dominance)
  colnames(variances) <- c("h2","d2","t2","q2")[1:dominance]
  j=1
  for (j in 1:dominance) {
    ans <- readBinMat(sub(pattern="X",replacement=c("a","d","t","q")[j],x="tmp/ETA_X_b.bin"))
    effect.mean[[j]] <- apply(ans,2,mean)
    effect.sd <- apply(ans,2,sd)
    Ne <- apply(ans,2,function(x){effectiveSize(mcmc(x))})
    effect.SE[[j]] <- effect.sd/sqrt(Ne)
    variances[,j] <- apply(tcrossprod(data@Z %*% data@geno[[k]][[j]],ans),2,var)
  }
  varE <- scan("tmp/varE.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]
  h2 <- variances/(apply(variances,1,sum) + varE)  #proportion of variance for each term
  Ne <- apply(h2,2,effectiveSize)
  return.var <- cbind(Estimate=round(apply(h2,2,mean),2),SE=round(sqrt(apply(h2,2,var)/Ne),2))
  
  effects <- vector("list",length=min(dominance,2))
  names(effects) <- c("add","dom")[1:min(dominance,2)]

  effects$add <- data.frame(Haplotype=attr(data@geno,"haplotypes"),Mean=effect.mean[[1]],SE=effect.SE[[1]],stringsAsFactors = F)
  effects$add <- effects$add[order(effects$add$Haplotype),]

  if (dominance > 1) {
    diplotypes <- strsplit(attr(data@geno,"diplotypes"),split="+",fixed=T)
    effects$dom <- data.frame(Haplotype1=sapply(diplotypes,function(x){x[1]}),Haplotype2=sapply(diplotypes,function(x){x[2]}),Mean=effect.mean[[2]],SE=effect.SE[[2]],stringsAsFactors = F)
    effects$dom <- effects$dom[order(effects$dom$Haplotype1,effects$dom$Haplotype2),]
  }

  #Plotting
  #Additive effects
  tmp <- strsplit(split=".",x=effects$add$Haplotype,fixed=T)
  founders <- sapply(tmp,function(x){x[1]})
  plotme <- data.frame(effects$add,founders=factor(founders))
  plotme$Haplotype <- rep(1:data@ploidy,length(levels(plotme$founders)))
  plotA <- ggplot(data=plotme,aes(x=Haplotype,y=Mean,fill=founders)) + 
    labs(title = paste("Trait:", trait),
         subtitle = paste("Marker:", marker)) + 
    theme_bw() + 
    scale_x_continuous(name="Haplotype",labels=1:data@ploidy,breaks=1:data@ploidy) +
    geom_bar(stat="identity",position="dodge",colour="black") + 
    ylab("Additive Effect") + 
    scale_fill_viridis_d(name="")+
    theme(text = element_text(size=13)) + 
    facet_grid(.~founders,scales = "free_x") +
    geom_errorbar(aes(ymax=Mean+1.96*SE, ymin = Mean-1.96*SE, width = 0.2))
  
  if (dominance > 1) {
    #Dominance effects
    plotme <- data.frame(x=c(effects$dom$Haplotype1,effects$dom$Haplotype2),y=c(effects$dom$Haplotype2,effects$dom$Haplotype1),z=rep(effects$dom$Mean,2)) #construct symmetric matrix
    plotme <- plotme[!duplicated(plotme[,1:2]),] #remove duplication of diagonals
    plotD <- ggplot(data=plotme,aes(x=x,y=y,fill=z)) +
      geom_tile() +
      scale_fill_viridis_c(name="Dominance\nEffect") +
      labs(title = paste("Trait:", trait),
         subtitle = paste("Marker:", marker)) +
      theme_bw() +
      xlab("") +
      ylab("") +
      theme(text = element_text(size=13),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
      coord_fixed(ratio=1)
      
    return(list(R2=ans1$R2,deltaDIC=ans1$DIC-ans0$DIC,var=return.var,effects=effects,plots=list(add=plotA,dom=plotD)))
  } else {
    return(list(R2=ans1$R2,deltaDIC=ans1$DIC-ans0$DIC,var=return.var,effects=effects,plots=list(add=plotA)))
  }
}

