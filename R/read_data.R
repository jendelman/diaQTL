#' Read data files
#' 
#' Reads genotype, pedigree, and phenotype data files 
#' 
#' First three columns of the genotype file are marker, chromosome, and position. Columns 4 through (n+4) correspond to the n individuals of the population. The genotype information for each marker x individual combination is a string with the format "state|state|state...=>prob|prob|prob...", where "state" refers to the genotype state and "prob" is the genotype probability in decimal format. Only states with nonzero probabilities need to be listed. The encoding for the states in tetraploids is described in the documentation for the F1codes and S1codes datasets that come with the package. For diploids, there are 4 F1 genotype codes, 1,2,3,4, which correspond to allele combinations 1-3,1-4,2-3,2-4, respectively; the S1 genotype codes 1,2,3 correspond to 1-1,1-2,2-2, respectively. For the phenotype file, first column is id, followed by traits, and then any fixed effects. Pass a character vector for the function argument "fixed" to specify whether each effect is a factor or numeric covariate. The number of traits is deduced based on the number of columns. Binary traits must be coded "N"/"Y" and are converted to 0/1 internally for analysis by probit regression. 
#'
#' @param genofile File with map and genotype probabilities 
#' @param ploidy Allowable values are 2 or 4
#' @param pedfile File with pedigree information (four column format: id,population,mother,father)
#' @param phenofile File with phenotype data (optional)
#' @param fixed If there are fixed effects, this is a character vector of "factor" or "numeric"
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
#' 
read_data <- function(genofile,ploidy=4,pedfile,phenofile=NULL,fixed=NULL) {
  
  dominance = TRUE
  data <- read.csv(genofile,as.is=T,check.names=F)
  gid <- colnames(data)[-(1:3)]
  dupe <- which(duplicated(gid))
  if (length(dupe) > 0) {
    stop("Duplicate names in genotype file.")
  }
  
  ped <- read.csv(pedfile,as.is=T)
  colnames(ped) <- c("id","population","mother","father")
  rownames(ped) <- ped[,1]
  missing <- setdiff(gid,ped$id)
  if(length(missing)>0) {
    stop("Not all genotyped individuals are in the pedigree file.")
  }
  n.gid <- length(gid)
  cat(paste("N =",n.gid,"individuals with pedigree and genotype information \n"))
  ped <- ped[which(ped$id %in% gid),]
  ped <- ped[match(gid,ped$id),]
  
  if (!is.null(phenofile)) {
    pheno <- read.csv(phenofile,check.names=F)
    pheno[,1] <- as.character(pheno[,1])
    missing2 <- setdiff(pheno[,1],gid)
    if(length(missing2)>0) {
      stop("Not all phenotyped individuals are in the genotype file")
    }
    
    n2 <- length(intersect(gid,pheno[,1]))
    cat(paste("N =",n2,"individuals with phenotype and genotype data \n"))
  }
  
  cat("Preparing genotype data...\n")
  genoX <- make_X(ped,ploidy=ploidy,dominance=dominance)
  
  map <- data[,1:3]
  chroms <- unique(map[,2])
  n.chrom <- length(chroms)
  first.chrom.marker <- match(chroms,map[,2])
  m <- nrow(map)
  n <- ncol(data)-3
  n.state <- ifelse(ploidy==2,4,100)
  
  genoA <- vector("list",length=m) #dim1 = id, dim2 = alleles
  names(genoA) <- map[,1]
  attr(genoA,"alleles") <- attr(genoX$A,"alleles")
  attr(genoA,"id") <- gid
  if (dominance) {
    genoD <- vector("list",length=m) #dim1 = id, dim2 = alleles
    names(genoD) <- map[,1]
    attr(genoD,"allele.pairs") <- attr(genoX$D,"allele.pairs")
    attr(genoD,"id") <- gid
  } else {
    genoD <- NULL
  }
  
  for (j in 1:m) {
    if (j %in% first.chrom.marker) {
      cat(paste("Chromosome",map[j,2],"\n"))
    }
    tmp <- unlist(lapply(data[j,3+1:n],strsplit,split="=>",fixed=T),recursive = F)
    tmp2 <- lapply(tmp,strsplit,split="|",fixed=T)
    states <- lapply(tmp2,function(x){as.integer(x[[1]])})
    ind <- rep(1:n,sapply(states,length))
    genoprob <- lapply(tmp2,function(x){y<-as.numeric(x[[2]])
    y/sum(y)})
    tmp3 <- mapply(FUN=function(u,v,w){
      tcrossprod(u[,v],t(w))
    },genoX$A,states,genoprob)
    if(class(tmp3)=="matrix"){ # in case of all probs=1
      tmp3 = split(tmp3, rep(1:ncol(tmp3), each = nrow(tmp3)))
    }
    genoA[[j]] <- Matrix(t(sapply(tmp3,function(x){as.vector(x)})),sparse=TRUE)
    
    if (dominance) {
      tmp3 <- mapply(FUN=function(u,v,w){
        tcrossprod(u[,v],t(w))
      },genoX$D,states,genoprob)
      if(class(tmp3)=="matrix"){ # in case of all probs=1
        tmp3 = split(tmp3, rep(1:ncol(tmp3), each = nrow(tmp3)))
      }
      genoD[[j]] <- Matrix(t(sapply(tmp3,function(x){as.vector(x)})),sparse=TRUE)
    }
  }
  
  data <- new(Class="diallel_geno",ploidy=as.integer(ploidy),ped=ped,map=map,geno=list(A=genoA,D=genoD))
  rm(list=c("genoA","genoD","ped","map"))
  tmp <- gc(verbose=FALSE)
  
  if (is.null(phenofile)) {
    return(data)
  } else {
    crossID <- data.frame(x=factor(data@ped$population))
    if(length(levels(crossID$x))==1){
      Xped <- sparse.model.matrix(~1,data=crossID)
    }else{
      Xped <- sparse.model.matrix(~x,data=crossID)
    }
    rownames(Xped) <- gid
    colnames(Xped) <- levels(crossID$x)
    
    pheno[,1] <- factor(pheno[,1],levels=gid)
    Z <- sparse.model.matrix(~id-1,data=pheno)
    colnames(Z) <- gid
    X <- Z%*%Xped
    
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
