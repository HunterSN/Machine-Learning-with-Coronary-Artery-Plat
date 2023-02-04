rm(list = ls())


symbol_index = read.csv('symbol_Index.csv',header = T, stringsAsFactors = F)
CAQ_index = read.csv('CAQ_index.csv',header = T, stringsAsFactors = F)
s_data = read.csv('source.csv',header = T, stringsAsFactors = F)
indexMM = as.data.frame(symbol_index$MM[which(symbol_index$MM != '无')])

sym = s_data[,c(1,34:87)]


symMM = as.data.frame(c())
catch = as.data.frame(indexMM)

for (i in 1:nrow(sym)) {
  print(i)
  catch[,2] = i
  catch2 = t(sym[i,c(2:55)])
  catch[,3] = catch2[,1]
  symMM = rbind(symMM,catch)
}
names(symMM) = c('index', 'num', 'type')

indexMM_CAQ = as.data.frame(indexMM)
names(indexMM_CAQ) = c('index')
indexMM_CAQ  = rbind(indexMM_CAQ,CAQ_index) 

sym_CAQ = s_data[,c(1,34:87,121:125)]

symMM_CAQ = as.data.frame(c())
catch = as.data.frame(indexMM_CAQ)

for (i in 1:nrow(sym_CAQ)) {
  print(i)
  catch[,2] = i
  catch2 = t(sym_CAQ[i,c(2:60)])
  catch[,3] = catch2[,1]
  symMM_CAQ = rbind(symMM_CAQ,catch)
}
names(symMM_CAQ) = c('index', 'num', 'type')

symMM = symMM[which(symMM$type != 0),]
symMM_CAQ = symMM_CAQ[which(symMM_CAQ$type != 0),]

symMM = symMM[,c(1,2)]
symMM_CAQ = symMM_CAQ[,c(1,2)]

#save.image(file = "sym&CAQ_Process.RData")


load('sym&CAQ_Process.RData')
library(openxlsx)
write.xlsx(symMM, file = "symMM.xlsx", colNames = TRUE)
write.xlsx(symMM_CAQ, file = "symMM_CAQ.xlsx", colNames = TRUE)


#save.image(file = "zgxfTapShaiApriori.RData")

#library(openxlsx)
#write.xlsx(shai, file = "zgxfZongApriori.xlsx", colNames = TRUE)

#load('knnData.RData')
#rm(list=ls()[which(ls()!='ZoutKNN')])



