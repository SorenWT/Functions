cvtest_matlab <- function(filein,fileout){

require('cvequality')
require('RJSONIO')
  
  library(cvequality)
  library(readxl)

  sheets <- readxl::excel_sheets(filein)
  filedata <- lapply(sheets, function(X) readxl::read_excel(filein, sheet = X))
  filedata <- lapply(filedata, as.data.frame)
  names(filedata) <- sheets
     
#filedata <- read.csv(filein)
  
  testoutput <- list()
  for(i in 1:length(filedata)){
    testoutput[[i]] <- with(filedata[[i]],asymptotic_test(X,G))
  }

testjson <- toJSON(testoutput)

write(testjson,file=fileout,sep="")

}