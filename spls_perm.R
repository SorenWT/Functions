spls_perm <- function(file1in,file2in,ncomps,nperms){
  
  require('RJSONIO')
  require('permuco')
  
  dat1 <- read.csv(file1in)
  dat1 = data.matrix(dat1)
  
  dat2 <- read.csv(file2in)
  dat2 = data.matrix(dat2)
  
  cvopts <- cv.spls(dat1,dat2,K=ncomps,eta=seq(0.1,0.9,0.1))
  
  plsmdl <- spls(dat1,dat2,K=cvopts$K.opt,eta=cvopts$eta.opt)
  
  for(i in 1:nperms){
    dat1perm = dat1[sample(1:nrows(dat1),nrows(dat1),replace=FALSE), ]
    cvperm <- cv.spls(dat1perm,dat2,K=ncomps,eta=seq(0.1,0.9,0.1))
    plsperm <- spls(dat1perm,dat2,K=cvperm$K.opt,eta=cvopts$eta.opt)
  }
  
  plsmdl
  
  mdljson <- toJSON(plsmdl)
  
  write(mdljson,file=paste(filein,'_results.json',sep=''))
}