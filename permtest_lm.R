permtest_lm <- function(filein,frmla){
  
  require('RJSONIO')
  require('permuco')

  dat <- read.csv(filein)
  
  mdl <- lmperm(frmla,data=dat,np=10000)
  
  mdljson <- toJSON(mdl)
  
  write(mdljson,file=paste(filein,'_results.json',sep=''))
  }