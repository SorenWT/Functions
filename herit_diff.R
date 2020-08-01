herit_diff <- function(mzin1,dzin1,mzin2,dzin2,mdltype){
  
  require('umx')
  #require('OpenMx')
  require('RJSONIO')
  require('XML')
  
  mzdat1 = read.csv(mzin1)
  dzdat1 = read.csv(dzin1)
  
  factnames = names(mzdat)
  factnames_short = sub('_T.*','',factnames)
  factnames_short = factnames_short[seq(from=1,to=length(factnames)-1,by=2)]
  
  mzdat = umx_scale_wide_twin_data(varsToScale=factnames_short,sep='_T',data=mzdat)
  dzdat = umx_scale_wide_twin_data(varsToScale=factnames_short,sep='_T',data=dzdat)
  
  mdl1 <- umxACE(selDVs=factnames_short,mzData=mzdat,dzData=dzdat,sep='_T',intervals=TRUE)
  
  stdmdl <- umxSummaryACE(mdl,returnStd=TRUE,digits=5)
  
  
  
  
  
  
}