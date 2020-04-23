artanova <- function(filein,frmla,fileout){
  require('RJSONIO')
  require('ARTool')
  require('emmeans')

  library('ARTool')  

  dat <- read.csv(filein)
  
  mdl <- art(frmla,data=dat)
  
  results <- anova(mdl)
  
  #results$part.eta.sq = with(results, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
  
  #library('emmeans')
  
  #multcompare <- list()
  
  # for(i in 1:nrow(results)){
  #   if(results$`Pr(>F)`[i] < 0.05){
  #     model.lm <- artlm(mdl,results$Term[i])
  #     marginal = emmeans(model.lm,results$Term[i])
  #     if(grepl(':',results$Term[i])){
  #       tmp <- contrast(marginal,method='pairwise',adjust='tukey')
  #       multcompare[[results$Term[i]]] <- summary(tmp)
  #     } else {
  # 
  #       tmp <- pairs(marginal,adjust='tukey')
  #       multcompare[[results$Term[i]]] <- summary(tmp)
  #     }    
  #   }
  # }
  
  resjson <- toJSON(results)
  
  #mcjson <- toJSON(multcompare)
  
  write(resjson,file=paste(fileout,'_results.json',sep=''))
  
  #write(mcjson,file=paste(fileout,'_multcompare.json',sep=''))
}
