#' @importFrom utils combn
#' 
make_X <- function(ped,ploidy,dominance) {
  # Xadd is used to map genotype probabilities for each F1 population to additive founder effects
  # ped is data frame with variables id,mother,father
  # Xadd has dimensions nind x (ploidy*nf) x (# F1 genotype states), where nf = number of founders and nind is number of individuals in the population
  # alleles have .1, .2, etc. appended to founder name
  # Xdom is used to map genotype probabilities for each population to digenic dominance founder effects
  # Xdom has dimensions nind x (# possible digenic effects) x (# F1 genotype states)
  
  founders <- as.character(unique(c(ped$mother,ped$father)))
  nf <- length(founders)
  nind <- nrow(ped)
  alleles <- as.character(1:(ploidy*nf))
  n.add <- length(alleles)
  
  digen <- apply(combn(alleles,2),2,paste,collapse="+")
  digen <- c(apply(cbind(alleles,alleles),1,paste,collapse="+"),digen)
  n.dom <- length(digen)
  n.state <- ifelse(ploidy==2,4,100)
  
  if (ploidy==2) {
    S1code <- rbind(c(1,1,2),c(1,2,2))
    F1code.mom <- c(1,1,2,2)
    F1code.dad <- c(3,4,3,4)
  } else {
    #tetraploid  
    tmp <- strsplit(S1codes$State,split="-",fixed=T)
    S1code <- sapply(tmp,function(x){as.integer(x)})
    tmp <- strsplit(F1codes$State,split="-",fixed=T)
    F1code.mom <- sapply(tmp,function(x){as.integer(x[1:2])})
    F1code.dad <- sapply(tmp,function(x){as.integer(x[3:4])})
  }
  
  Xadd <- vector("list",length=nind) #dim1 = alleles, dim2 = genotype states
  if (dominance) {
    Xdom <- vector("list",length=nind) #dim1 = allele.pairs, dim2 = genotype states
  } else {
    Xdom <- NULL 
  }
  
  mom <- match(ped$mother,founders)
  dad <- match(ped$father,founders)
  for (i in 1:nind) {
    km <- (mom[i]-1)*ploidy
    kd <- (dad[i]-1)*ploidy
    if (km==kd) { 
      #self
      genocode <- S1code
      for (j in 1:ploidy) {
        genocode[genocode==j] <- km+j
      }
    } else {
      #F1 pop
      mom.code <- F1code.mom
      dad.code <- F1code.dad
      for (j in 1:ploidy) {
        mom.code[mom.code==j] <- km+j
        dad.code[dad.code==(j+ploidy)] <- kd+j
      }
      genocode <- rbind(mom.code,dad.code)
    }
    x <- split(genocode,col(genocode))  #matrix to vector coercion is columnwise
    tmp <- lapply(x,table)
    names(tmp) <- NULL
    tmp2 <- unlist(tmp)
    ix <- match(names(tmp2),alleles)
    Xadd[[i]] <- sparseMatrix(i=ix,j=rep(1:length(tmp),times=sapply(tmp,length)),x=as.integer(tmp2),dims=c(n.add,n.state))
    if (dominance) {
      genocode2 <- apply(genocode,2,function(x){
          y <- combn(x,2)
          apply(y,2,function(z){paste(sort(z),collapse="+")})
        })
      if (ploidy==2) {
        genocode2 <- matrix(genocode2,nrow=1)
      }
      x <- split(genocode2,col(genocode2))  #matrix to vector coercion is columnwise
      tmp <- lapply(x,table)
      names(tmp) <- NULL
      tmp2 <- unlist(tmp)
      ix <- match(names(tmp2),digen)
      Xdom[[i]] <- sparseMatrix(i=ix,j=rep(1:length(tmp),times=sapply(tmp,length)),x=as.integer(tmp2),dims=c(n.dom,n.state))
    }
  }
  
  tmp <- expand.grid(1:ploidy,founders,stringsAsFactors = F)
  alleles <- apply(cbind(tmp[,2],tmp[,1]),1,paste,collapse=".")
  attr(Xadd,"alleles") <- alleles
  
  if (dominance) {
    digen <- apply(combn(alleles,2),2,paste,collapse="+")
    digen <- c(apply(cbind(alleles,alleles),1,paste,collapse="+"),digen)
    attr(Xdom,"allele.pairs") <- digen
  }
  return(list(A=Xadd,D=Xdom))
}
