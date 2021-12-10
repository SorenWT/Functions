ranova_matlab <- function(filein,DV,WID,withinvar,fileout)
  
library('rstatix')
library('RJSONIO')
library('tidyverse')

mydata <- read.csv(filein)
mydata <- as_tibble(mydata)

res.aov <- anova_test(mydata, dv=DV,wid=WID,within=withinvar)

anovatbl <- get_anova_table(res.aov)

anovajson <- toJSON(res.aov) 
anovatbljson <- toJSON(anovatbl)

write(anovajson,file=paste(fileout,'_mdl.json',sep=''))
write(anovatbljson,file=paste(fileout,'_tbl.json',sep=''))




