#' Read data files
#' 
#' Reads genotype, pedigree, and phenotype data files 
#' 
#' The first 3 columns of the genotype file should be the genetic map (labeled marker, chrom, cM), and a fourth column for a reference genome position (labeled bp) can also be included. The map is followed by the members of the population. The genotype data for each marker x individual combination is a string with the format "state|state|state...=>prob|prob|prob...", where "state" refers to the genotype state and "prob" is the genotype probability in decimal format. Only states with nonzero probabilities need to be listed. The encoding for the states in tetraploids is described in the documentation for the F1codes and S1codes datasets that come with the package. For diploids, there are 4 F1 genotype codes, 1,2,3,4, which correspond to haplotype combinations 1-3,1-4,2-3,2-4, respectively; the S1 genotype codes 1,2,3 correspond to 1-1,1-2,2-2, respectively. For the phenotype file, first column is id, followed by traits, and then any fixed effects. Pass a character vector for the function argument "fixed" to specify whether each effect is a factor or numeric covariate. The number of traits is deduced based on the number of columns. Binary traits must be coded N/Y and are converted to 0/1 internally for analysis by probit regression. Missing data in the phenotype file should be coded as NA. The parameter \code{dominance} specifies the maximum value of dominance that can be used in subsequent analysis: 1 = additive, 2 = digenic dominance, 3 = trigenic dominance, 4 = quadrigenic dominance. The default is dominance = ploidy, which allows the full range of dominance models in functions such as \code{\link{scan1}} and \code{\link{fitQTL}}, but this requires the most RAM. Output files from the BGLR package are stored in a folder named 'tmp' in the current directory.
#'
#' @param genofile File with map and genotype probabilities 
#' @param ploidy Either 2 or 4
#' @param pedfile File with pedigree data (id,parent1,parent2)
#' @param phenofile File with phenotype data (optional)
#' @param fixed If there are fixed effects, this is a character vector of "factor" or "numeric"
#' @param bin.markers TRUE/FALSE whether to bin markers with the same cM position
#' @param dominance Maximum value of dominance that will be used for analysis. Default = ploidy.
#' @param n.core Number of cores for parallel execution
#' 
#' @return Variable of class \code{\link{diallel_geno}} if phenofile is NULL, otherwise \code{\link{diallel_geno_pheno}}
#' 
#' @examples
#' \dontrun{
#'   ## Get the location of raw csv files examples
#'   genocsv = system.file( "vignette_data", "potato_geno.csv", package = "diaQTL" )
#'   pedcsv = system.file( "vignette_data", "potato_ped.csv", package = "diaQTL" )
#'   phenocsv = system.file( "vignette_data", "potato_pheno.csv", package = "diaQTL" )
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
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport
#' 
read_data <- function(genofile,ploidy=4,pedfile,phenofile=NULL,
                      fixed=NULL,bin.markers=TRUE,dominance=NULL,n.core=1) {
  stopifnot(ploidy %in% c(2,4))
  if (is.null(dominance)) {
    dominance <- ploidy
  }
  if ((dominance > 2) & (ploidy==2)) {
    stop("Only digenic dominance exists for diploids.")
  }

  data <- read.csv(genofile,as.is=T,check.names=F)
  cM <- grep("cM",colnames(data),fixed=T)
  bp <- grep("bp",colnames(data),fixed=T)
  if (length(cM)==0) {
    stop("Map requires column labeled cM")
  }
  if (length(bp)>0) {
    map <- data[,c(1:2,cM,bp)]
    colnames(map) <- c("marker","chrom","cM","bp")
  } else {
    map <- data[,c(1:2,cM)]
    colnames(map) <- c("marker","chrom","cM")
  }
  
  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)
  m <- nrow(map)
  if (bin.markers) {
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
  
  id <- colnames(data)[ncol(map):ncol(data)]
  dupe <- which(duplicated(id))
  if (length(dupe) > 0) {
    stop("Duplicate names in genotype file.")
  }
  
  ped <- read.csv(pedfile,as.is=T)
  if (ncol(ped) > 3) {
    stop("Pedigree file should have 3 columns: id,parent1,parent2")
  }
  colnames(ped) <- c("id","parent1","parent2")
  rownames(ped) <- ped[,1]
  id <- intersect(id,ped$id)
  cat(paste(length(id),"individuals with pedigree and genotype data \n"))

  if (!is.null(phenofile)) {
    pheno <- read.csv(phenofile,check.names=F,stringsAsFactors = T)
    pheno[,1] <- as.character(pheno[,1])
    pheno <- pheno[pheno[,1] %in% id,]
    id <- intersect(id,pheno[,1])
    cat(paste(length(id),"individuals with pedigree, genotype and phenotype data \n"))
  } 
  
  ped <- ped[match(id,ped$id),]
  data <- as.matrix(data[,match(id,colnames(data))])
  
  #GCA
  parents <- sort(unique(c(ped$parent1,ped$parent2)))
  tmp <- data.frame(parent1=factor(ped$parent1,levels=parents,ordered=T),parent2=factor(ped$parent2,levels=parents,ordered=T))
  if (length(parents)==1) {
    X.GCA <- Matrix(2,nrow=nrow(tmp),ncol=1)
    colnames(X.GCA) <- parents
  } else {
    X.GCA <- sparse.model.matrix(~parent1-1,tmp) + sparse.model.matrix(~parent2-1,tmp)
    colnames(X.GCA) <- gsub(pattern="parent1",replacement="",colnames(X.GCA))
  }
  rownames(X.GCA) <- ped$id
  
  cat("Preparing genotype data...\n")
  genoX <- make_X(ped,ploidy,dominance)
  n.state <- ifelse(ploidy==2,4,100)
  
  f1 <- function(j,data,genoX,ploidy,dominance) {
    tmp <- unlist(lapply(data[j,],strsplit,split="=>",fixed=T),recursive = F)
    tmp2 <- lapply(tmp,strsplit,split="|",fixed=T)
    states <- lapply(tmp2,function(x){as.integer(x[[1]])})
    genoprob <- lapply(tmp2,function(x){y<-as.numeric(x[[2]])
        y/sum(y)})
    
    geno <- vector("list",length=dominance)
    for (q in 1:dominance) {
      tmp3 <- mapply(FUN=function(u,v,w){tcrossprod(u[,v],t(w))},u=genoX[[q]],v=states,w=genoprob)
      if(is.list(tmp3)){ 
        tmp3 <- sapply(tmp3,function(x){as.vector(x)})
      } 
      geno[[q]] <- Matrix(t(tmp3))
    }
    return(geno)
  }
  f2 <- function(ix,geno,ploidy) {
    id <- attr(geno,"id")
    n <- length(id)
    K <- matrix(0,nrow=n,ncol=n)
    dimnames(K) <- list(id,id)
    m <- length(ix)
    for (i in ix) {
      K <- K + as.matrix(tcrossprod(geno[[i]][[1]]))
    }
    return(K/m/ploidy)
  }
  
  cl <- makeCluster(n.core)
  clusterExport(cl=cl,varlist=NULL)
  
  geno <- parLapply(cl, bin.ix, f1, data=data,genoX=genoX,ploidy=ploidy,dominance=dominance)
  #geno <- lapply(bin.ix,f1, data=data,genoX=genoX,ploidy=ploidy,dominance=dominance)
  names(geno) <- bin.names
  attr(geno,"id") <- id
  attr(geno,"haplotypes") <- attr(genoX,"haplotypes")
  if (dominance > 1) {
    attr(geno,"diplotypes") <- attr(genoX,"diplotypes")
  }
  
  ### additive relationship matrix for polygenic effects
  map2 <- map[map$marker %in% bin.names,]
  markers <- split(map2$marker,f=map2$chrom)
  A <- parLapply(cl,markers,f2,geno=geno,ploidy=ploidy)
  names(A) <- chroms
  
  stopCluster(cl)

  data <- new(Class="diallel_geno",ploidy=as.integer(ploidy),input=data[bin.ix,],Xa=genoX[[1]],
              dominance=as.integer(dominance),X.GCA=X.GCA,map=map,geno=geno,A=A)
  rownames(data@input) <- bin.names
  
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
