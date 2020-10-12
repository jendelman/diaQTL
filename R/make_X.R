#' @importFrom arrangements combinations
#' 
make_X <- function(ped,ploidy,dominance) {
  # ped is data frame with variables id,parent1,parent2
  # X1 is used to map genotype probabilities for each F1 population to additive founder effects
  # X1 has dimensions nind x (ploidy*nf) x (# F1 genotype states), where nf = number of founders and nind is number of individuals in the population
  # X2 is used to map genotype probabilities for each population to digenic dominance founder effects
  # X2 has dimensions nind x (# possible digenic effects) x (# F1 genotype states)
  # X3 and X4 are for trigenic and quadrigenic effects
  # haplotypes have .1, .2, etc. appended to founder name
  
  founders <- sort(unique(c(ped$parent1,ped$parent2)))
  crosses <- unique(apply(ped[,2:3],1,function(x){x <- sort(match(x,founders))
                                                  paste(x,collapse="x")}))
  n.cross <- length(crosses)
  nf <- length(founders)
  nind <- nrow(ped)
  haplotypes <- as.character(1:(ploidy*nf))
  n1 <- length(haplotypes)
  
  X <- vector("list",length=dominance)
  for (q in 1:dominance) {
    X[[q]] <- vector("list",length=nind)
  }
  
  if (dominance > 1) {
    #digenic
    gen2 <- character(0)
    
    if (ploidy==4) {
      #uniparental terms
      for (p1 in 1:nf) {
        gen2 <- c(gen2,apply(combinations(x=ploidy*(p1-1)+1:ploidy,k=2,replace=T),1,paste,collapse="+"))
      }
    }
    
    #biparental terms
    for (i in 1:n.cross) {
      p1 <- as.integer(substr(crosses[i],1,1))
      p2 <- as.integer(substr(crosses[i],3,3))
      if ((p1!=p2) | (ploidy==2)) {
        gen2 <- c(gen2,apply(expand.grid(x=ploidy*(p1-1)+1:ploidy,y=ploidy*(p2-1)+1:ploidy),1,paste,collapse="+"))
      }
    }
    n2 <- length(gen2)
  }

  if (dominance > 2) {
    #trigenic = digenic from one parent plus monogenic from the other
    gen3 <- character(0)
    for (i in 1:n.cross) {
      p1 <- as.integer(substr(crosses[i],1,1))
      p2 <- as.integer(substr(crosses[i],3,3))
      
      p1.di <- combinations(x=ploidy*(p1-1)+1:ploidy,k=2,replace=T)
      tmp <- as.matrix(expand.grid(x=1:10,y=ploidy*(p2-1)+1:ploidy))
      gen3 <- c(gen3,apply(tmp,1,function(z){paste(sort(c(p1.di[z[1],],z[2])),collapse="+")}))

      if (p1!=p2) {
        p2.di <- combinations(x=ploidy*(p2-1)+1:ploidy,k=2,replace=T)
        tmp <- as.matrix(expand.grid(x=1:10,y=ploidy*(p1-1)+1:ploidy))
        gen3 <- c(gen3,apply(tmp,1,function(z){paste(sort(c(p2.di[z[1],],z[2])),collapse="+")}))
      }
    }
    n3 <- length(gen3)
  }
  
  if (dominance > 3) {
    #quadrigenic = digenic from both parents
    gen4 <- character(0)
    for (i in 1:n.cross) {
      p1 <- as.integer(substr(crosses[i],1,1))
      p2 <- as.integer(substr(crosses[i],3,3))
    
      p1.di <- combinations(x=ploidy*(p1-1)+1:ploidy,k=2,replace=T)
      p2.di <- combinations(x=ploidy*(p2-1)+1:ploidy,k=2,replace=T)
      tmp <- as.matrix(expand.grid(x=1:10,y=1:10))
      gen4 <- c(gen4,apply(tmp,1,function(z){paste(sort(c(p1.di[z[1],],p2.di[z[2],])),collapse="+")}))
    }
    n4 <- length(gen4)
  }
  
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

  mom <- match(ped$parent1,founders)
  dad <- match(ped$parent2,founders)
  i=1
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
    ix <- match(names(tmp2),haplotypes)
    X[[1]][[i]] <- sparseMatrix(i=ix,j=rep(1:length(tmp),times=sapply(tmp,length)),x=as.integer(tmp2),dims=c(n1,n.state))
    
    if (dominance > 1) {
      #digenic
      genocode2 <- apply(genocode,2,function(x){
        y <- combinations(x,2,replace=F)
        apply(y,1,function(z){paste(sort(z),collapse="+")})
        })
      if (ploidy==2) {
        genocode2 <- matrix(genocode2,nrow=1)
      } 
      x <- split(genocode2,col(genocode2))  #matrix to vector coercion is columnwise
      tmp <- lapply(x,table)
      names(tmp) <- NULL
      tmp2 <- unlist(tmp)
      ix <- match(names(tmp2),gen2)
      X[[2]][[i]] <- sparseMatrix(i=ix,j=rep(1:length(tmp),times=sapply(tmp,length)),x=as.integer(tmp2),dims=c(n2,n.state))
    }
    
    if (dominance > 2) {
      #trigenic
      genocode3 <- apply(genocode,2,function(x){
        y <- combinations(x,3,replace=F)
        apply(y,1,function(z){paste(sort(z),collapse="+")})
      })
      x <- split(genocode3,col(genocode3))  #matrix to vector coercion is columnwise
      tmp <- lapply(x,table)
      names(tmp) <- NULL
      tmp2 <- unlist(tmp)
      ix <- match(names(tmp2),gen3)
      X[[3]][[i]] <- sparseMatrix(i=ix,j=rep(1:length(tmp),times=sapply(tmp,length)),x=as.integer(tmp2),dims=c(n3,n.state))
    }
    
    if (dominance > 3) {
      #quadrigenic
      genocode4 <- apply(genocode,2,function(x){
        paste(sort(x),collapse="+")
      })
      ix <- match(genocode4,gen4)
      X[[4]][[i]] <- sparseMatrix(i=ix,j=1:length(ix),x=rep(1,length(ix)),dims=c(n4,n.state))
    }
  }
  
  tmp <- expand.grid(1:ploidy,founders,stringsAsFactors = F)
  haplotypes <- apply(cbind(tmp[,2],tmp[,1]),1,paste,collapse=".")
  attr(X,"haplotypes") <- haplotypes
  
  if (dominance > 1) {
    
    gen2 <- character(0)
    if (ploidy==4) {
      #uniparental terms
      for (p1 in 1:nf) {
        gen2 <- c(gen2,apply(combinations(x=haplotypes[ploidy*(p1-1)+1:ploidy],k=2,replace=T),1,paste,collapse="+"))
      }
    }
    
    #biparental terms
    for (i in 1:n.cross) {
      p1 <- as.integer(substr(crosses[i],1,1))
      p2 <- as.integer(substr(crosses[i],3,3))
      if ((p1!=p2) | (ploidy==2)) {
        gen2 <- c(gen2,apply(expand.grid(x=haplotypes[ploidy*(p1-1)+1:ploidy],y=haplotypes[ploidy*(p2-1)+1:ploidy]),1,paste,collapse="+"))
      }
    }
    attr(X,"diplotypes") <- gen2
  }

  return(X)
}
