#' Fit multiple QTL model
#' 
#' Fit multiple QTL model
#' 
#' Argument \code{qtl} is a data frame with columns `marker` and `dominance` 
#' to specify the marker name and highest order effect (1 = additive, 2 = digenic dominance, 
#' 3 = trigenic dominance, 4 = quadrigenic dominance). All effects up to the value in `dominance` 
#' are included.  Optional argument \code{epistasis} is a data frame with columns `marker1` and 
#' `marker2`, where each row specifies an additive x additive epistatic interaction. 
#' The number of burn-in and total iterations in \code{params} can be estimated using 
#' \code{\link{set_params}}. Parameter \code{CI.prob} sets the probability (e.g., 0.90, 0.95) 
#' for the Bayesian credible interval for the estimated effects (to disable plotting of the CI, 
#' use \code{CI.prob=NULL}). 
#' 
#' @param data variable of class \code{\link{diallel_geno_pheno}}
#' @param trait name of trait
#' @param qtl data frame, see Details
#' @param epistasis optional data frame, see Details
#' @param polygenic TRUE/FALSE whether to include additive polygenic effect
#' @param params list containing the number of burn-in (burnIn) and total iterations (nIter)
#' @param CI.prob probability for Bayesian credible interval
#' 
#' @return List containing
#' \describe{
#' \item{deltaDIC}{DIC relative to model with GCA but no QTL effects}
#' \item{resid}{residuals}
#' \item{var}{matrix with proportion of variance for the effects}
#' \item{effects}{list with two matrices, `additive` and `digenic`, with markers on the rows and effects on the columns}
#' \item{plots}{list of ggplot objects, one for each marker, containing elements `additive` and `digenic`. The digenic plot has digenic effects above the diagonal and the sum of additive and digenic effects below the diagonal.}
#' }
#' @examples
#' \dontrun{
#' ## getting minimum burnIn and nIter for one qtl
#' set_params(data = diallel_example, 
#'            trait = "tuber_shape", 
#'            q = 0.05, 
#'            r = 0.025, 
#'            qtl = data.frame(marker="solcap_snp_c2_25522",dominance=2),
#'            polygenic = TRUE)
#'            
#' ## additive effects
#' fit1 <- fitQTL(data = diallel_example, 
#'                trait = "tuber_shape", 
#'                params = list(burnIn=100,nIter=5000), 
#'                qtl = data.frame(marker="solcap_snp_c2_25522",dominance=1),
#'                CI.prob = 0.9)
#'
#' ## additive + digenic dominance effects                            
#' fit2 <- fitQTL(data = diallel_example, 
#'                trait = "tuber_shape", 
#'                params = list(burnIn=100,nIter=5000), 
#'                qtl = data.frame(marker="solcap_snp_c2_25522",dominance=2),
#'                CI.prob=0.9)
#'                
#' ## getting minimum burnIn and nIter for two qtl with epistasis
#' set_params(data = diallel_example, 
#'            trait = "tuber_shape", 
#'            q = 0.05, 
#'            r = 0.025, 
#'            qtl = data.frame(marker=c("PotVar0099535","solcap_snp_c2_25522"),
#'                             dominance=c(2,1)),
#'            epistasis = data.frame(marker1="solcap_snp_c2_25522",marker2="PotVar0099535"),
#'            polygenic = TRUE)
#'            
#' ## additive + digenic dominance effects for both QTL
#' fit3 <- fitQTL(data = diallel_example, trait = "tuber_shape", 
#'                params = list(burnIn=100,nIter=5000),
#'                qtl = data.frame(marker=c("PotVar0099535","solcap_snp_c2_25522"),
#'                                 dominance=c(2,2)), 
#'                polygenic = TRUE, CI.prob = 0.9)
#'                
#' ## additive + digenic dominance effects for both QTL + their epistatic effects
#' fit4 <- fitQTL(data = diallel_example, trait = "tuber_shape", 
#'                params = list(burnIn=100,nIter=5000),
#'                qtl = data.frame(marker=c("PotVar0099535","solcap_snp_c2_25522"),
#'                                 dominance=c(2,2)), 
#'                epistasis = data.frame(marker1="solcap_snp_c2_25522",marker2="PotVar0099535"),
#'                polygenic = TRUE, CI.prob = 0.9)
#'                
#' ## additive + digenic dominance effects for three QTL + all their epistatic effects
#' fit5 <- fitQTL(data = diallel_example, trait = "tuber_shape", 
#'                params = list(burnIn=100,nIter=5000),
#'                qtl = data.frame(marker=c("PotVar0099535",
#'                                          "solcap_snp_c1_6427",
#'                                          "solcap_snp_c2_25522"),
#'                                 dominance=c(2,2,2)), 
#'                epistasis = data.frame(marker1=c("solcap_snp_c2_25522",
#'                                                 "solcap_snp_c2_25522",
#'                                                 "PotVar0099535"),
#'                                       marker2=c("PotVar0099535",
#'                                                 "solcap_snp_c1_6427",
#'                                                 "solcap_snp_c1_6427")),
#'                polygenic = TRUE, CI.prob = 0.9)
#'
#' }
#'                  
#' @export
#' 
#' @importFrom stats var sd quantile
#' @import ggplot2
#' @importFrom BGLR readBinMat
#' @importFrom rlang .data

fitQTL <- function(data,trait,qtl,epistasis=NULL,polygenic=FALSE,params=list(burnIn=100,nIter=5000),CI.prob=0.9) {
  
  stopifnot(inherits(data,"diallel_geno_pheno"))
  stopifnot(trait %in% colnames(data@pheno))
  stopifnot(qtl$marker %in% data@map$marker)
  if (max(qtl$dominance) > data@dominance) {
    stop("Dominance degree cannot exceed value used with read_data.")
  }
  
  markers <- get_bin(qtl$marker,data@map)
  n.qtl <- length(markers)
  Xqtl <- vector("list",n.qtl)
  for (i in 1:n.qtl) {
    dominance <- qtl$dominance[i]
    Xqtl[[i]] <- vector("list",dominance)
    for (j in 1:dominance) {
      Xqtl[[i]][[j]] <- data@geno[[markers[i]]][[j]]  #data@Z %*% 
    }
  }
  
  if (!is.null(epistasis)) {
    marker1 <- get_bin(epistasis$marker1,data@map)
    marker2 <- get_bin(epistasis$marker2,data@map)
    if (length(setdiff(c(marker1,marker2),markers)) > 0) {
      stop("Markers in epistasis should also be in qtl")
    }
    n.epi <- nrow(epistasis)
    Xaa <- vector("list",n.epi)
    for (i in 1:n.epi) {
      Xaa[[i]] <- faa(j=marker1[i],k=marker2[i],data=data)  #data@Z %*% 
    }
  } else {
    n.epi <- 0
    marker1 <- marker2 <- Xaa <- NULL
  }

  if (polygenic) {
    poly.marker <- unique(c(markers,marker1,marker2))
    chroms <- setdiff(unique(data@map$chrom),data@map$chrom[match(poly.marker,data@map$marker)])
    if (length(chroms)==0) {
      stop("There are no chromosomes remaining for the polygenic effect.")
    }
    polyG <- IBDmat(data=data,dominance=0,chrom=chroms) 
  } else {
    polyG <- NULL
  }
  
  y <- data@pheno[,trait]
  if (inherits(y,"factor")) {
    response <- "ordinal"
  } else {
    response <- "gaussian"
  }
  
  if (params$burnIn < 0) {
    #undocumented option, used by set_params
    set.params <- TRUE
    params$burnIn <- 0
  } else {
    set.params <- FALSE
  }
  params <- list(response=response,nIter=params$nIter,burnIn=params$burnIn)

  if (!is.null(CI.prob) && CI.prob < 0) {
    #undocumented option, used by scan1
    ans1 <- runBGLR(params=params,y=y,Xfix=data@X,Z=data@Z,Xqtl=Xqtl,Xaa=Xaa,saveEffects=FALSE)
    return(ans1)
  } else {
    #null hypothesis
    ans0 <- runBGLR(params=params,y=y,Xfix=data@X,Z=data@Z,Xgca=data@X.GCA,saveEffects=FALSE)
    #alternate hypothesis
    ans1 <- runBGLR(params=params,y=y,Xfix=data@X,Z=data@Z,Xqtl=Xqtl,Xaa=Xaa,polyG=polyG,saveEffects=TRUE)
  }
  deltaDIC <- ans1$DIC - ans0$DIC
  
  trace.length <- params$nIter - params$burnIn
  ans <- get_trace(qtl=qtl,epistasis=epistasis,polygenic=polygenic,params=params)
  variances <- NULL
  haplotypes <- attr(data@geno,"haplotypes")
  n.hap <- length(haplotypes)
  additive.effects <- matrix(as.numeric(NA),nrow=n.qtl,ncol=n.hap)
  colnames(additive.effects) <- haplotypes
  rownames(additive.effects) <- qtl$marker
  CI.lower <- CI.upper <- additive.effects
  
  if (data@dominance > 1) {
    diplotypes <- attr(data@geno,"diplotypes")
    n.diplo <- length(diplotypes)
    digenic.effects <- matrix(as.numeric(NA),nrow=n.qtl,ncol=n.diplo)
    colnames(digenic.effects) <- diplotypes
    rownames(digenic.effects) <- qtl$marker
  }
    
  for (i in 1:n.qtl) {
    dominance <- qtl$dominance[i]
    tmp <- matrix(0,nrow=trace.length,ncol=dominance)
    colnames(tmp) <- paste(qtl$marker[i],c("additive","digenic","trigenic","quadrigenic")[1:dominance])
    for (j in 1:dominance) {
      tmp[,j] <- apply(tcrossprod(Xqtl[[i]][[j]],ans$qtl[[i]][[j]]),2,var)
    }
    variances <- cbind(variances,tmp)
    
    ## additive and digenic effects
    additive.effects[i,] <- apply(ans$qtl[[i]][[1]],2,mean)
    if (!is.null(CI.prob)) {
      tmp2 <- apply(ans$qtl[[i]][[1]],2,quantile,p=c(0.5-CI.prob/2,0.5+CI.prob/2))
      CI.lower[i,] <- tmp2[1,]
      CI.upper[i,] <- tmp2[2,]
    }
    if (dominance > 1) {
      digenic.effects[i,] <- apply(ans$qtl[[i]][[2]],2,mean)
    }
  }
  
  if (data@dominance > 1) {
    diplo2 <- strsplit(diplotypes,split="+",fixed=T)
    diplo2 <- data.frame(hap1=sapply(diplo2,function(x){x[1]}),hap2=sapply(diplo2,function(x){x[2]}))
    
    #remove digenic effects between copies of the same haplotype unless that parent was selfed
    max.dosage <- apply(data@X.GCA,2,max)
    selfed <- setdiff(names(max.dosage),names(which(max.dosage==1)))
    parent1 <- sapply(strsplit(diplo2[,1],split=".",fixed=T),"[",1)
    keep <- which(diplo2[,1]!=diplo2[,2] | parent1 %in% selfed)
    digenic.effects <- matrix(digenic.effects[,keep],nrow=n.qtl)
    colnames(digenic.effects) <- diplotypes[keep]
    rownames(digenic.effects) <- qtl$marker
    diplo2 <- diplo2[keep,]
  }
  
  if (!is.null(epistasis)) {
    n.epi <- nrow(epistasis)
    tmp <- matrix(0,nrow=trace.length,ncol=n.epi)
    colnames(tmp) <- paste(paste(epistasis$marker1,epistasis$marker2,sep="+"),"epistasis")
    epistasis.effects <- matrix(NA,nrow=ncol(ans$epi[[1]]),ncol=n.epi)
    for (i in 1:n.epi) {
      tmp[,i] <- apply(tcrossprod(Xaa[[i]],ans$epi[[i]]),2,var)
      epistasis.effects[,i] <- apply(ans$epi[[i]],2,mean)
    }
    colnames(epistasis.effects) = paste(paste(epistasis$marker1,epistasis$marker2,sep="+"))
    rownames(epistasis.effects) = apply(expand.grid(attr(data@geno,"haplotypes"),
                                                    attr(data@geno,"haplotypes"))[,2:1],1,
                                        paste0,collapse="+")
    variances <- cbind(variances,tmp)
  } 
  
  if (polygenic) {
    variances <- cbind(variances,polygenic=ans$poly*mean(diag(polyG)))
  }  
  if (params$response=="gaussian") {
    variances <- cbind(variances,residual=ans$resid)
    h2 <- variances/apply(variances,1,sum)
  } else {
    y2 <- as.integer(y)-1
    R2 <- 1 - sum(ans1$resid^2,na.rm=T)/sum((y2-mean(y2,na.rm=T))^2,na.rm=T)
    h2 <- variances/apply(variances,1,sum)*R2
  }
  
  if (set.params) {
    return(variances)
  }
  
  if (!is.null(CI.prob) & response=="gaussian") {
    tmp <- apply(h2,2,quantile,p=c(0.5-CI.prob/2,0.5+CI.prob/2))
    return.var <- cbind(Mean=round(apply(h2,2,mean),2),CI.lower=round(tmp[1,],2),CI.upper=round(tmp[2,],2))
  } else {
    return.var <- cbind(Mean=round(apply(h2,2,mean),2))
  }
  
  #Now make plots
  plots <- vector("list",n.qtl)
  names(plots) <- qtl$marker
  for (i in 1:n.qtl) {
    dominance <- qtl$dominance[i]
    plots[[i]] <- vector("list",length=min(dominance,2))
    names(plots[[i]]) <- c("additive","digenic")[1:min(dominance,2)]
  
    tmp <- strsplit(split=".",x=haplotypes,fixed=T)
    parents <- sapply(tmp,"[",1)
    hap.num <- as.integer(sapply(tmp,"[",2))
    plot.data <- data.frame(parent=factor(parents),x=hap.num,
                            mean=additive.effects[i,],
                            CI.lower=CI.lower[i,],CI.upper=CI.upper[i,])
    plots[[i]]$additive <- ggplot(data=plot.data,aes(x=.data$x,y=.data$mean)) + 
      labs(title = paste("Trait:", trait),subtitle = paste("Marker:", qtl$marker[i])) + 
      theme_bw() + 
      scale_x_continuous(name="Haplotype",labels=1:data@ploidy,breaks=1:data@ploidy) +
      geom_bar(stat="identity",position="dodge",colour="black",fill="grey50") + 
      ylab("Additive Effect") + theme(text = element_text(size=13)) + 
      facet_grid(.~parent,scales = "free_x")
    
    if (!is.null(CI.prob)) {
      plots[[i]]$additive <- plots[[i]]$additive + geom_errorbar(aes(ymax=.data$CI.upper, ymin=.data$CI.lower, width = 0.2))
    }
    
    if (dominance > 1) {
      plot.data <- data.frame(x=c(diplo2$hap1,diplo2$hap2),y=c(diplo2$hap2,diplo2$hap1),
                              z=c(digenic.effects[i,],digenic.effects[i,]+additive.effects[i,diplo2$hap1]+additive.effects[i,diplo2$hap2]))
      plot.data <- plot.data[!duplicated(plot.data[,1:2]),]
      plots[[i]]$digenic <- ggplot(data=plot.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + 
        geom_tile() + scale_fill_gradient2(name="") + 
        labs(title = paste("Trait:", trait),
         subtitle = paste("Marker:", qtl$marker[i])) +
        theme_bw() + xlab("") + ylab("") +
        theme(text = element_text(size=13),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
        coord_fixed(ratio=1)
    }
  }
  
  if (max(qtl$dominance)>1) {
    effects <- list(additive=t(additive.effects),digenic=t(digenic.effects))
  } else {
    effects <- list(additive=t(additive.effects))
  }
  
  if(!is.null(epistasis)){
    effects[['epistasis']] = epistasis.effects
    for(i in 1:n.epi){
      n.hap = length(attr(data@geno,"haplotypes"))
      plot.data <- data.frame(x=rep(attr(data@geno,"haplotypes"),each=n.hap),
                              y=rep(attr(data@geno,"haplotypes"),n.hap),
                              z=c(effects$epistasis[,i]))
      plots$epistasis[[i]] <- ggplot(data=plot.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + 
        geom_tile() + scale_fill_gradient2(name="") + 
        labs(title = paste("Trait:", trait),
             subtitle = "Epistasis") +
        theme_bw() + xlab(epistasis[i,1]) + ylab(epistasis[i,2]) +
        theme(text = element_text(size=13),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
        coord_fixed(ratio=1)
    }
    names(plots$epistasis) = colnames(effects$epistasis)
  }
  return(list(deltaDIC=deltaDIC,resid=ans1$resid,
              var=return.var,effects=effects,plots=plots))
}

