#' Fit a single QTL model
#' 
#' Fit a single QTL model
#' 
#' Standard errors of the posterior mean estimates are calculated by dividing the SD of the Markov Chain by the square root of the effective number of iterations, which is calculated by function \code{effectiveSize} in R package \code{coda}. The error bars on the plot of additive effects correspond to +/- 1.96*SE (95 percent confidence interval). For binary traits, R2 = the squared phi correlation. The additive and dominance variances are reported as a proportion of the total variance: h2=Va/(Va+Vd+Vresid) and d2=Vd/(Va+Vd+Vresid). Parameter \code{dominance} controls the genetic model: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. 
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
#' \item{h2}{Mean and SE for proportion of variance due to additive effects}
#' \item{effectsA}{Mean and SE of the additive effects for parental haplotypes}
#' \item{plotA}{ggplot object for additive effects}
#' }
#' If dominance > 1, the list also contains
#' \describe{
#' \item{d2}{Mean and SE for proportion of variance due to dominance (all orders)}
#' \item{effectsD}{Mean and SE of the digenic dominance effects}
#' \item{plotD}{ggplot object for digenic dominance effects}
#' }
#' 
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
  
  add_b <- readBinMat('ETA_add_b.bin')
  haplotype.mean <- apply(add_b,2,mean)
  haplotype.SD <- apply(add_b,2,sd)
  haplotype.Ne <- apply(add_b,2,function(x){effectiveSize(mcmc(x))})
  effectsA <- data.frame(Haplotype=attr(data@geno,"haplotypes"),Mean=haplotype.mean,SE=haplotype.SD/sqrt(haplotype.Ne),stringsAsFactors = F)
  effectsA <- effectsA[order(effectsA$Haplotype),]

  #Make plot of additive effects
  tmp <- strsplit(split=".",x=effectsA$Haplotype,fixed=T)
  founders <- sapply(tmp,function(x){x[1]})
  plotme <- data.frame(effectsA,founders=factor(founders))
  plotme$Haplotype <- 1:data@ploidy
  plotA <- ggplot(data=plotme,aes(x=Haplotype,y=Mean,fill=founders)) + 
    labs(title = paste("Trait:", trait),
         subtitle = paste("Marker:", marker)) + 
    theme_bw() + 
    xlab("Haplotype") + 
    geom_bar(stat="identity",position="dodge",colour="black") + 
    ylab("Additive Effect") + 
    scale_fill_viridis_d(name="")+
    theme(text = element_text(size=13),
          axis.ticks.x=element_blank()) + 
    facet_grid(.~founders,scales = "free_x") +
    geom_errorbar(aes(ymax=Mean+1.96*SE, ymin = Mean-1.96*SE, width = 0.2))

  if (dominance > 1) {
    dom_b <- readBinMat('ETA_dom_b.bin')
    digenic.mean <- apply(dom_b,2,mean)
    digenic.SD <- apply(dom_b,2,sd)
    digenic.Ne <- apply(dom_b,2,function(x){effectiveSize(mcmc(x))})
    HaplotypePair <- strsplit(attr(data@geno,"haplotype.pairs"),split="+",fixed=T)
    effectsD <- data.frame(Haplotype1=sapply(HaplotypePair,function(x){x[1]}),Haplotype2=sapply(HaplotypePair,function(x){x[2]}),Mean=digenic.mean,SE=digenic.SD/sqrt(digenic.Ne),stringsAsFactors = F)
    effectsD <- effectsD[order(effectsD$Haplotype1,effectsD$Haplotype2),]

    #plot of dominance effects
    plotme <- data.frame(x=c(effectsD$Haplotype1,effectsD$Haplotype2),y=c(effectsD$Haplotype2,effectsD$Haplotype1),z=rep(effectsD$Mean,2)) #construct symmetric matrix
    plotme <- plotme[-which(duplicated(plotme[,1:2])),] #remove duplication of diagonals
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
  }
  
  #Heritability
  varA <- apply(tcrossprod(data@Z %*% data@geno[[k]][[1]],add_b),2,var)
  varE <- scan("varE.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]
  
  if (dominance==1) {
    h2 <- varA/(varA+varE)
    h2 <- c(mean(h2),as.numeric(sqrt(var(h2)/effectiveSize(h2))))
    names(h2) <- c("Estimate","SE")
    return(list(R2=ans1$R2,deltaDIC=ans1$DIC-ans0$DIC,h2=h2,effectsA=effectsA,plotA=plotA))
    
  } else {
    varD <- apply(tcrossprod(data@Z %*% data@geno[[k]][[2]],dom_b),2,var)
    if (dominance > 2) {
      d3_b <- readBinMat('ETA_d3_b.bin')
      varD <- varD + apply(tcrossprod(data@Z %*% data@geno[[k]][[3]],d3_b),2,var)
    }
    if (dominance > 3) {
      d4_b <- readBinMat('ETA_d4_b.bin')
      varD <- varD + apply(tcrossprod(data@Z %*% data@geno[[k]][[4]],d4_b),2,var)
    }
    h2 <- varA/(varA+varD+varE)
    h2 <- c(mean(h2),as.numeric(sqrt(var(h2)/effectiveSize(h2))))
    names(h2) <- c("Estimate","SE")

    d2 <- varD/(varA+varD+varE)
    d2 <- c(mean(d2),as.numeric(sqrt(var(d2)/effectiveSize(d2))))
    names(d2) <- c("Estimate","SE")

    return(list(R2=ans1$R2,deltaDIC=ans1$DIC-ans0$DIC,h2=h2,d2=d2,effectsA=effectsA,effectsD=effectsD,plotA=plotA,plotD=plotD))
  }
}

