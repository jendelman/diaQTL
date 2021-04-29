#' Generate diaQTL input files from RABBIT MagicReconstruct 
#' 
#' Generate diaQTL input files from RABBIT MagicReconstruct
#' 
#' @param rabbit.outfile name of RABBIT output file
#' @param ped.file name of RABBIT pedigree file
#' @param outstem prefix for the pedigree and genotype files for diaQTL
#' 
#' @return NULL
#' @export
#' @importFrom utils write.csv

convert_rabbit <- function(rabbit.outfile,ped.file,outstem) {
  con <- file(ped.file,"r",)
  suppressWarnings(temp <- readLines(con))
  close(con)
  ix <- grep("Pedigree-Information",temp)
  ped <- temp[(ix[1]+2):(ix[2]-1)]
  tmp2 <- strsplit(ped,split=",",fixed=T)
  gen <- as.integer(sapply(tmp2,"[",1))
  stopifnot(max(gen)==1)  #only F1 or S1 pops allowed
  
  pop <- as.integer(sapply(tmp2,"[",2))
  parent1 <- as.integer(sapply(tmp2,"[",4))
  parent2 <- as.integer(sapply(tmp2,"[",5))
  
  tmp2 <- strsplit(temp[(ix[2]+2):length(temp)],split=",",fixed=T)
  progeny.id <- sapply(tmp2,"[",1)
  pop.id <- as.integer(sapply(tmp2,"[",2))
  
  ped <- data.frame(id=progeny.id,
                    parent1=parent1[match(pop.id,pop)],
                    parent2=parent2[match(pop.id,pop)])
  
  con <- file(rabbit.outfile,"r",)
  suppressWarnings(temp <- readLines(con))
  close(con)
  
  tmp2 <- strsplit(temp[2:4],split=",",fixed=T)
  map <- data.frame(marker=tmp2[[1]][-1],
                    chrom=tmp2[[2]][-1],
                    cM=round(as.numeric(tmp2[[3]][-1]),2))
  
  ix <- grep("magicReconstruct",temp)
  tmp2 <- temp[5:(ix[2]-1)]
  tmp2 <- strsplit(tmp2,split=",",fixed=T)
  haplotypes <- sapply(tmp2,"[",1)
  parents <- unique(gsub("_Paternal","",gsub("_Maternal","",haplotypes)))
  ped$parent1 <- parents[ped$parent1]
  ped$parent2 <- parents[ped$parent2]
  write.csv(ped,file=paste0(outstem,"diaQTL_ped.csv"),row.names=F)
  
  #genotype file
  k <- min(grep("genotype",temp[ix],ignore.case=F))
  x <- strsplit(temp[c((ix[k]+2):(ix[k+1]-1))],split=",")
  ng <- length(x)
  geno.code <- sapply(x,"[",2)
  uni.pop <- unique(pop.id)
  n.pop <- length(uni.pop)
  
  #for each population, construct mapping from geno.code to diaQTL codes
  genotypes <- vector("list",n.pop)
  names(genotypes) <- uni.pop
  for (i in 1:n.pop) {
    j <- match(uni.pop[i],pop)
    tmp2 <- c(parent1[j],parent2[j])
    p1 <- max(tmp2)
    p2 <- setdiff(tmp2,p1)
    p1.1 <- 2*(p1-1)+1
    p1.2 <- p1.1 + 1
    if (length(p2)>0) {
      p2.1 <- 2*(p2-1)+1
      p2.2 <- p2.1 + 1
      tmp2 <- c(paste(p1.1,p2.1,sep="|"),paste(p1.1,p2.2,sep="|"),
                paste(p1.2,p2.1,sep="|"),paste(p1.2,p2.2,sep="|"))
    } else {
      tmp2 <- c(paste(p1.1,p1.1,sep="|"),paste(p1.1,p1.2,sep="|"),paste(p1.2,p1.2,sep="|"))
    }
    genotypes[[i]] <- match(geno.code,tmp2)
  }

  k <- grep("genoprob",temp[ix],ignore.case=F)
  x <- strsplit(temp[c((ix[k]+4):(ix[k+1]-1))],split=",")
  y <- sapply(x,function(x){round(as.numeric(x[-1]),3)})
  m <- nrow(map)
  n <- ncol(y)/ng
  w <- strsplit(sapply(x,"[",1),split="_genotype",fixed=T)
  out1 <- matrix(0,nrow=m,ncol=n)
  colnames(out1) <- unique(sapply(w,"[",1))
  for (i in 1:m) {
    z <- split(y[i,],rep(1:n,each=ng))
    for (j in 1:n) {
      geno.code <- genotypes[[as.character(pop.id[j])]]
      ix <- which(z[[j]] > 0)
      stopifnot(!is.na(geno.code[ix]))
      out1[i,j] <- paste(paste(geno.code[ix],collapse="|"),
                         paste(z[[j]][ix],collapse="|"),sep="=>")
    }
  }
  write.csv(cbind(map,out1),file=paste0(outstem,"diaQTL_geno.csv"),row.names=F)
}