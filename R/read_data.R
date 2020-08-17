#' Read data files
#' 
#' Reads genotype, pedigree, and phenotype data files 
#' 
#' The first 3 or 4 columns of the genotype file are the map (marker, chrom, bp and/or cM), followed by the members of the population. The genotype information for each marker x individual combination is a string with the format "state|state|state...=>prob|prob|prob...", where "state" refers to the genotype state and "prob" is the genotype probability in decimal format. Only states with nonzero probabilities need to be listed. The encoding for the states in tetraploids is described in the documentation for the F1codes and S1codes datasets that come with the package. For diploids, there are 4 F1 genotype codes, 1,2,3,4, which correspond to haplotype combinations 1-3,1-4,2-3,2-4, respectively; the S1 genotype codes 1,2,3 correspond to 1-1,1-2,2-2, respectively. For the phenotype file, first column is id, followed by traits, and then any fixed effects. Pass a character vector for the function argument "fixed" to specify whether each effect is a factor or numeric covariate. The number of traits is deduced based on the number of columns. Binary traits must be coded N/Y and are converted to 0/1 internally for analysis by probit regression. Parameter \code{dominance} controls the genetic model: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance.
#'
#' @param genofile File with map and genotype probabilities 
#' @param ploidy Either 2 or 4
#' @param pedfile File with pedigree information (id,mother,father)
#' @param phenofile File with phenotype data (optional)
#' @param fixed If there are fixed effects, this is a character vector of "factor" or "numeric"
#' @param bin.markers TRUE/FALSE whether to bin markers with the same cM position
#' @param dominance Dominance degree (1-4). See Details.
#' @param n.core Number of cores for parallel execution (only available from Linux or Mac command line)

#' 
#' @return Variable of class \code{\link{diallel_geno}} if phenofile is NULL, otherwise \code{\link{diallel_geno_pheno}}
#' 
#' @examples
#' \dontrun{
#'   ## Get the location of raw csv files examples
#'   genocsv = system.file( "tutorial", "potato_geno.csv", package = "diaQTL" )
#'   pedcsv = system.file( "tutorial", "potato_ped.csv", package = "diaQTL" )
#'   phenocsv = system.file( "tutorial", "potato_pheno.csv", package = "diaQTL" )
#'   
#'   ## Check their location in the system
#'   print(genocsv)
#'   print(pedcsv)
#'   print(phenocsv)
#'   
#'   ## Load them in R
#'   diallel_example <- read_data(genofile = genocsv,
#'                                ploidy = 4,
#'                                pedfile = pedcsv,
#'                                phenofile = phenocsv)
#' }
#' 
#' @export
#' @import Matrix
#' @importFrom utils read.csv
#' @importFrom methods new
#' @importFrom parallel mclapply
#' 
read_data <- function(genofile,ploidy=4,pedfile,phenofile=NULL,fixed=NULL,bin.markers=T,dominance=2,n.core=1) {
  
  if ((dominance > 2) & (ploidy==2)) {
    stop("Only digenic dominance exists for diploids.")
  }
  
  data <- read.csv(genofile,as.is=T,check.names=F)
  cM <- grep("cM",colnames(data))
  bp <- grep("bp",colnames(data))
  if ((length(cM)>0) & (length(bp)>0)) {
    map <- data[,c(1:2,cM,bp)]
  }
  if ((length(cM)>0) & (length(bp)==0)) {
    map <- data[,c(1:2,cM)]
  }
  if ((length(cM)==0) & (length(bp)>0)) {
    map <- data[,c(1:2,bp)]
  }
  if ((length(cM)==0) & (length(bp)==0)) {
    stop("Invalid map position. Must be cM and/or bp.")
  }
  colnames(map) <- c("marker","chrom",colnames(map)[-(1:2)])
  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)
  m <- nrow(map)
  if ((length(cM)>0) & bin.markers) {
    tmp <- apply(map[,c("chrom","cM")],1,paste,collapse="_")
    bins <- unique(tmp)
    map$bin <- match(tmp,bins)
    n.bin <- length(bins)
    bin.ix <- match(1:n.bin,map$bin)
    bin.names <- map$marker[bin.ix]
  } else {
    map$bin <- 1:m
    n.bin <- m
    bin.ix <- 1:m
    bin.names <- map$marker[bin.ix]
  }
  first.chrom.marker <- map$marker[match(chroms,map$chrom)]
  
  id <- colnames(data)[-(1:max(cM,bp))]
  dupe <- which(duplicated(id))
  if (length(dupe) > 0) {
    stop("Duplicate names in genotype file.")
  }
  
  ped <- read.csv(pedfile,as.is=T)
  colnames(ped) <- c("id","mother","father")
  rownames(ped) <- ped[,1]
  missing <- setdiff(id,ped$id)
  if(length(missing)>0) {
    stop("Not all genotyped individuals are in the pedigree file.")
  }
  cat(paste("N =",length(id),"individuals with pedigree and genotype information \n"))
  ped <- ped[which(ped$id %in% id),]
  ped <- ped[match(id,ped$id),]
  
  #GCA
  parents <- sort(unique(c(ped$mother,ped$father)))
  tmp <- data.frame(mother=factor(ped$mother,levels=parents,ordered=T),father=factor(ped$father,levels=parents,ordered=T))
  if (length(parents)==1) {
    X.GCA <- Matrix(1,nrow=nrow(tmp),ncol=1)
    colnames(X.GCA) <- parents
  } else {
    X.GCA <- sparse.model.matrix(~mother-1,tmp)/2 + sparse.model.matrix(~father-1,tmp)/2
    colnames(X.GCA) <- gsub(pattern="mother",replacement="",colnames(X.GCA))
  }
  rownames(X.GCA) <- ped$id
  
  
  if (!is.null(phenofile)) {
    pheno <- read.csv(phenofile,check.names=F)
    pheno[,1] <- as.character(pheno[,1])
    missing2 <- setdiff(pheno[,1],id)
    if(length(missing2)>0) {
      stop("Not all phenotyped individuals are in the genotype file")
    }
    n2 <- length(intersect(id,pheno[,1]))
    cat(paste("N =",n2,"individuals with phenotype and genotype data \n"))
  }
  
  cat("Preparing genotype data...\n")
  genoX <- make_X(ped,ploidy,dominance)
  n.state <- ifelse(ploidy==2,4,100)
  
  f1 <- function(j,data,genoX,id,ploidy,dominance) {
    tmp <- unlist(lapply(data[j,id],strsplit,split="=>",fixed=T),recursive = F)
    tmp2 <- lapply(tmp,strsplit,split="|",fixed=T)
    states <- lapply(tmp2,function(x){as.integer(x[[1]])})
    genoprob <- lapply(tmp2,function(x){y<-as.numeric(x[[2]])
        y/sum(y)})
    
    geno <- vector("list",length=dominance)
    for (q in 1:dominance) {
      tmp3 <- mapply(FUN=function(u,v,w){tcrossprod(u[,v],t(w))},u=genoX[[q]],v=states,w=genoprob)
      if(class(tmp3)=="matrix"){ # in case of all probs=1
        tmp3 = split(tmp3, rep(1:ncol(tmp3), each = nrow(tmp3)))
      }
      geno[[q]] <- Matrix(t(sapply(tmp3,function(x){as.vector(x)})),sparse=TRUE)
    }
    return(geno)
  }
  
  geno <- mclapply(X=bin.ix,FUN=f1,data=data,genoX=genoX,id=id,ploidy=ploidy,dominance=dominance,mc.cores=n.core)
  
  names(geno) <- bin.names
  attr(geno,"id") <- id
  attr(geno,"haplotypes") <- attr(genoX,"haplotypes")
  if (dominance > 1) {
    attr(geno,"haplotype.pairs") <- attr(genoX,"haplotype.pairs")
  }

  data <- new(Class="diallel_geno",ploidy=as.integer(ploidy),dominance=as.integer(dominance),X.GCA=X.GCA,map=map,geno=geno)
  
  if (is.null(phenofile)) {
    return(data)
  } else {
    pheno[,1] <- factor(pheno[,1],levels=id)
    Z <- sparse.model.matrix(~id-1,data=pheno)
    colnames(Z) <- id
    X <- sparse.model.matrix(~1,data=pheno)
    
    colp <- colnames(pheno)
    n.fixed <- length(fixed)
    n.trait <- length(colp)-n.fixed-1
    
    for (j in 1:n.trait) {
      if (inherits(pheno[,1+j],"factor")) {
        v <- levels(pheno[,1+j])
        if (length(v)!=2 | v[1]!="N" | v[2]!="Y") {
          stop(paste(colp[j+1],"detected as binary trait but levels are not N/Y"))
        }
      }
    }
    
    if (n.fixed > 0) {
      for (j in 1:n.fixed) {
        if (fixed[j]=="factor") {
          pheno[,1+n.trait+j] <- factor(pheno[,1+n.trait+j]) 
        }
      }
      fixed.vars <- colp[1+n.trait+1:n.fixed]
      a <- paste("~",paste(fixed.vars,collapse="+"),sep="")
      X2 <- sparse.model.matrix(eval(parse(text=a)),pheno)
      X <- cbind(X,X2[,-1])  #remove intercept
    }
    return(new(Class="diallel_geno_pheno",data,pheno=pheno[,1+0:n.trait],X=X,Z=Z))
  }
}
