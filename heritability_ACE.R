heritability_ACE <- function(mzin,dzin,fileout){
  
  require('umx')
  #require('OpenMx')
  require('RJSONIO')
  require('XML')
  
  mzdat = read.csv(mzin)
  dzdat = read.csv(dzin)
  
  factnames = names(mzdat)
  factnames_short = sub('_T.*','',factnames)
  factnames_short = factnames_short[seq(from=1,to=length(factnames)-1,by=2)]
  
  mzdat = umx_scale_wide_twin_data(varsToScale=factnames_short,sep='_T',data=mzdat)
  dzdat = umx_scale_wide_twin_data(varsToScale=factnames_short,sep='_T',data=dzdat)
  
  mdl <- umxACE(selDVs=factnames_short,mzData=mzdat,dzData=dzdat,sep='_T',intervals=TRUE)
  
  stdmdl <- umxSummaryACE(mdl,returnStd=TRUE,digits=5)
  
  acesmry <-list()
  acesmry$a <- stdmdl@submodels$top@matrices$a@values
  acesmry$c <- stdmdl@submodels$top@matrices$c@values
  acesmry$e <- stdmdl@submodels$top@matrices$e@values
  
  ci <- confint(stdmdl,run=TRUE)
  acesmry$ci <- ci
  
  smryfull <- summary(stdmdl)
  
  
  # Model reduction
  mdl_ce <- umxModify(stdmdl,update='a_r1c1',name='CE')
  stdmdl_ce <- umxSummaryACE(mdl_ce,returnStd=TRUE,digits=5)
  mdl_ae <- umxModify(stdmdl,update='c_r1c1',name='AE')
  stdmdl_ae <- umxSummaryACE(mdl_ae,returnStd=TRUE,digits=5)
  
  mdlcompare <- umxCompare(stdmdl,c(stdmdl_ae,stdmdl_ce),compareWeightedAIC=TRUE)
  
  if(which.min(mdlcompare$AIC)==1){
    mdl_reduc <- stdmdl
    }else if(which.min(mdlcompare$AIC)==2){
      mdl_reduc <- stdmdl_ae
    }else if(which.min(mdlcompare$AIC)==3){
      mdl_reduc <- stdmdl_ce
    }
  
  #mdl_reduc <- umxReduce(stdmdl,report='html')
  
  #mdl_reduc <- umxSummaryACE(mdl_reduc,returnStd=TRUE,digits=5)
  #mdlcompare <- readHTMLTable('tmp.html')
  
  reducsmry <-list()
  reducsmry$a <- mdl_reduc@submodels$top@matrices$a@values
  reducsmry$c <- mdl_reduc@submodels$top@matrices$c@values
  reducsmry$e <- mdl_reduc@submodels$top@matrices$e@values
  
  ci <- confint(mdl_reduc,run=TRUE)
  reducsmry$ci <- ci
  
  smryfull_reduc <- summary(stdmdl)
  
  
  # convert to JSON and write
  acesmryjson <- toJSON(acesmry)
  smryfulljson <- toJSON(smryfull)
  
  reducsmryjson <- toJSON(reducsmry)
  smryfull_reducjson <- toJSON(smryfull_reduc)
  
  mdlcomparejson <- toJSON(mdlcompare)
  
  write(acesmryjson,file=paste(fileout,'_acesmry.json',sep=''))
  write(smryfulljson,file=paste(fileout,'_fullsmry.json',sep=''))
  
  write(reducsmryjson,file=paste(fileout,'_reducsmry.json',sep=''))
  write(smryfull_reducjson,file=paste(fileout,'_fullsmry_reduc.json',sep=''))
  
  write(mdlcomparejson,file=paste(fileout,'_mdlcompare.json',sep=''))
}