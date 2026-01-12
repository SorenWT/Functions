sem_matlab <- function(filein,fileout,syntax,clustervar){
    
  library('lavaan')
  library('RJSONIO')
  
  print(syntax)

  mydata <- read.csv(filein)
  
  if(clustervar=='none'){ 
    thisfit <- sem(syntax,data=mydata)
  }else{
    thisfit <- sem(syntax,data=mydata,cluster=clustervar)
  }
  #thisfit <- sem(syntax,data=mydata,se="bootstrap")
  thissummary <- summary(thisfit)
  #paramscov <- vcov(thisfit)
  
  fitjson <- toJSON(thisfit) 
  summaryjson <- toJSON(thissummary)
  #covjson <- toJSON(params)
  
  write(fitjson,file=paste(fileout,'_mdl.json',sep=''))
  write(summaryjson,file=paste(fileout,'_summary.json',sep=''))
  #write(covjson,file=paste(fileout,'_cov.json',sep=''))
  }
