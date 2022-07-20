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
  haplotypes <- gsub("Paternal","p",gsub("Maternal","m",haplotypes))

  ped$parent1 <- parents[ped$parent1]
  ped$parent2 <- parents[ped$parent2]
  write.csv(ped,file=paste(outstem,"diaQTL_ped.csv",sep="_"),row.names=F)
  
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
  i=1
  F1code <- expand.grid(3:4,1:2)[,c(2,1)]
  S1code <- rbind(c(1,1),c(1,2),c(2,2))
    
  for (i in 1:n.pop) {
    j <- match(uni.pop[i],pop)
    p1 <- parents[parent1[j]]
    p2 <- parents[parent2[j]]
    
    if (p1!=p2) {
      q <- match(c(paste(p1,c("m","p"),sep="_"),paste(p2,c("m","p"),sep="_")),
                 haplotypes)
      tmp2 <- cbind(q[F1code[,1]],q[F1code[,2]])
    } else {
      q <- match(paste(p1,c("m","p"),sep="_"),haplotypes)
      tmp2 <- cbind(q[S1code[,1]],q[S1code[,2]])
    }
    tmp3 <- apply(tmp2,1,function(z){paste(sort(z,decreasing=T),collapse="|")})
    genotypes[[i]] <- match(geno.code,tmp3)
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
  write.csv(cbind(map,out1),file=paste(outstem,"diaQTL_geno.csv",sep="_"),row.names=F)
}
