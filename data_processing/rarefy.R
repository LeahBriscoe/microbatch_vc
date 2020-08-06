rarefy <- function(x,maxdepth){
  
  # x = totalsource
  # maxdepth = COVERAGE
  
  if(is.null(maxdepth)) return(x)
  
  if(!is.element(class(x), c('matrix', 'data.frame','array')))
    x <- matrix(x,nrow=nrow(x))
  nr <- nrow(x)
  nc <- ncol(x)
  
  for(i in 1:nrow(x)){
    if(sum(x[i,]) > maxdepth){
      prev.warn <- options()$warn
      options(warn=-1)
      s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
      options(warn=prev.warn)
      x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
    }
  }
  return(x)
}


test = t(otu_table[50:52,1:4])
rarefy(test,10)


