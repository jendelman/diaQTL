get_bin <- function(marker,map) {
  k1 <- match(marker,map$marker,nomatch=0)
  apply(array(k1),1,function(k) {
    if (k > 0) {
      return(map$marker[map$bin==map$bin[k]][1])
    } else {
      return(NA)
    }
  })
}
