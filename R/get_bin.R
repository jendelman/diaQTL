get_bin <- function(marker,map) {
  k <- match(marker,map$marker,nomatch=0)
  if (k > 0) {
    return(map$marker[map$bin==map$bin[k]][1])
  } else {
    return(NA)
  }
}
