diptest_matlab <- function(filein,fileout){
  
  require(diptest)
  require(RJSONIO)
  
  library(diptest)
  library(RJSONIO)
  
  dat <- read.csv(filein)
  
  dt <- lapply(dat,dip.test)
  
  testjson <- toJSON(dt)
  
  write(testjson,file=fileout,sep="")
  
  
}
  
  