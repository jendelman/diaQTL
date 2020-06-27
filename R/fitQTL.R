#' Fit a single QTL model
#' 
#' Fit a single QTL model
#' 
#' Standard errors of the posterior mean estimates are calculated by dividing the SD of the Markov Chain by the square root of the effective number of iterations, which is calculated by function \code{effectiveSize} in R package \code{coda}. The error bars on the plot of additive effects correspond to +/- 1.96*SE (95 percent confidence interval). For binary traits, R2 = the squared phi correlation. The additive and digenic dominance variances are reported as a proportion of the total variance: h2=Va/(Va+Vd+Vresid) and d2=Vd/(Va+Vd+Vresid). 
#' 
#' @param data Variable of class \code{\link{diallel_geno_pheno}}
#' @param trait Name of trait
#' @param marker Name of marker to fit as QTL
#' @param cofactor Name of marker to fit as cofactor
#' @param params List containing the number of burn-in (burnIn) and total iterations (nIter)
#' @param dominance Logical variable whether to include digenic dominance effects
#' 
#' @return List containing
#' \describe{
#' \item{R2}{Coefficient of determination}
#' \item{DIC}{Deviance Information Criterion}
#' \item{h2}{Mean and SE for proportion of variance due to additive effects}
#' \item{effectsA}{Mean and SE of the additive effects for parental alleles}
#' \item{plotA}{ggplot object for additive effects}
#' }
#' If dominance=T the list also contains
#' \describe{
#' \item{d2}{Mean and SE for proportion of variance due to digenic dominance effects}
#' \item{effectsD}{Mean and SE of the dominance effects for parental allele pairs}
#' \item{plotD}{ggplot object for dominance effects}
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
#'                  dominance = TRUE)
#'                  
#' }
#'                  
#' @export
#' 
#' @importFrom stats var sd
#' @importFrom coda effectiveSize
#' @import ggplot2
#' @importFrom BGLR readBinMat

fitQTL <- function(data,trait,marker,params,dominance=FALSE,cofactor=NULL) {
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  stopifnot(marker %in% names(data@geno$A))
  
  if (is.null(data@geno$D) & dominance) {
    stop("Dominance was FALSE in read_data")
  }

  y <- data@pheno[,trait]
  if (inherits(y,"factor")) {
    response <- "ordinal"
  } else {
    response <- "gaussian"
  }
  params <- list(response=response,nIter=params$nIter,burnIn=params$burnIn)

  if (!is.null(cofactor)) {
    Xcof <- data@Z%*%data@geno$A[[cofactor]]
  } else {
    Xcof <- NULL
  }
  k <- match(marker,data@map[,1])
  chrom <- as.character(data@map[k,2])

  #no marker model
  ans0 <- qtl1(y=y,X=data@X,Z=data@Z,params=params,Xcof=Xcof)
  
  #with marker
  if (!dominance) {
    ans1 <- qtl1(y=y,X=data@X,Z=data@Z,genoA=data@geno$A[[k]],params=params,Xcof=Xcof)
  } else {
    ans1 <- qtl1(y=y,X=data@X,Z=data@Z,genoA=data@geno$A[[k]],genoD=data@geno$D[[k]],params=params,Xcof=Xcof)
  }
  
  add_b <- readBinMat('ETA_add_b.bin')
  allele.mean <- apply(add_b,2,mean)
  allele.SD <- apply(add_b,2,sd)
  allele.Ne <- apply(add_b,2,function(x){effectiveSize(mcmc(x))})
  effectsA <- data.frame(Allele=attr(data@geno$A,"alleles"),Mean=allele.mean,SE=allele.SD/sqrt(allele.Ne),stringsAsFactors = F)
  effectsA <- effectsA[order(effectsA$Allele),]

  #Make plot of additive effects
  tmp <- strsplit(split=".",x=effectsA$Allele,fixed=T)
  founders <- sapply(tmp,function(x){x[1]})
  plotme <- data.frame(effectsA,founders=factor(founders))
  plotme$Allele <- 1:data@ploidy
  plotA <- ggplot(data=plotme,aes(x=Allele,y=Mean,fill=founders)) + 
    labs(title = paste("Trait:", trait),
         subtitle = paste("Marker:", marker)) + 
    theme_bw() + 
    xlab("Parental Allele") + 
    geom_bar(stat="identity",position="dodge",colour="black") + 
    ylab("Additive Effect") + 
#    scale_fill_brewer(name="",palette="Blues") + 
    scale_fill_viridis_d(name="")+
    theme(text = element_text(size=13),
          axis.ticks.x=element_blank()) + 
    facet_grid(.~founders,scales = "free_x") +
    geom_errorbar(aes(ymax=Mean+1.96*SE, ymin = Mean-1.96*SE, width = 0.2))

  if (dominance) {
    dom_b <- readBinMat('ETA_dom_b.bin')
    digenic.mean <- apply(dom_b,2,mean)
    digenic.SD <- apply(dom_b,2,sd)
    digenic.Ne <- apply(dom_b,2,function(x){effectiveSize(mcmc(x))})
    AllelePair <- strsplit(attr(data@geno$D,"allele.pairs"),split="+",fixed=T)
    effectsD <- data.frame(Allele1=sapply(AllelePair,function(x){x[1]}),Allele2=sapply(AllelePair,function(x){x[2]}),Mean=digenic.mean,SE=digenic.SD/sqrt(digenic.Ne),stringsAsFactors = F)
    effectsD <- effectsD[order(effectsD$Allele1,effectsD$Allele2),]
    
    #plot of dominance effects
    plotme <- data.frame(x=c(effectsD$Allele1,effectsD$Allele2),y=c(effectsD$Allele2,effectsD$Allele1),z=rep(effectsD$Mean,2)) #construct symmetric matrix
    plotme <- plotme[-which(duplicated(plotme[,1:2])),] #remove duplication of diagonals
    plotD <- ggplot(data=plotme,aes(x=x,y=y,fill=z)) + 
      geom_tile() + 
      scale_fill_viridis_c(name="Dominance\nEffect") + 
      labs(title = paste("Trait:", trait),
           subtitle = paste("Marker:", marker)) + 
      theme_bw() + 
      xlab("") + 
      ylab("") + 
      theme(text = element_text(size=13),axis.text.x = element_text(angle = 90)) + 
      coord_fixed(ratio=1)
  }
  
  #Heritability
  varA <- apply(tcrossprod(data@Z %*% data@geno$A[[k]],add_b),2,var)
  varE <- scan("varE.dat",quiet = T)[params$burnIn+1:(params$nIter-params$burnIn)]
  
  if (!dominance) {
    h2 <- varA/(varA+varE)
    h2 <- list(Mean=mean(h2),SE=as.numeric(sqrt(var(h2)/effectiveSize(h2))))
    return(list(R2=ans1$R2,deltaDIC=ans1$DIC-ans0$DIC,h2=h2,effectsA=effectsA,plotA=plotA))
  } else {
    varD <- apply(tcrossprod(data@Z %*% data@geno$D[[k]],dom_b),2,var)
    h2 <- varA/(varA+varD+varE)
    h2 <- list(Mean=mean(h2),SE=as.numeric(sqrt(var(h2)/effectiveSize(h2))))
  
    d2 <- varD/(varA+varD+varE) 
    d2 <- list(Mean=mean(d2),SE=as.numeric(sqrt(var(d2)/effectiveSize(d2))))
    return(list(R2=ans1$R2,deltaDIC=ans1$DIC-ans0$DIC,h2=h2,d2=d2,effectsA=effectsA,effectsD=effectsD,plotA=plotA,plotD=plotD))
  }
}

