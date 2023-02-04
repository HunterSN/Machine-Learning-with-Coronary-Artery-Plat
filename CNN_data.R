rm(list = ls())

load('sym&CAQ_Process.RData')

#load('RWR_PPIN_MM_target_CAQ_seedCAQ.RData')


#names(indexMM) = c('labels')
#names(CAQ_index) = c('labels')

#indexMM_Out = merge(indexMM, RWR_PPIN_MM_target_CAQ_seedCAQ, all.x = T)
#CAQ_index_Out = merge(CAQ_index, RWR_PPIN_MM_target_CAQ_seedCAQ, all.x = T)
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

CAQ = s_data[,c(1,121:125)]

symMM_CNN = merge(symMM, CAQ, all.x = T)

symMM_CNN = symMM_CNN[,-8]


MM_CAQ_catch = as.data.frame(c())
catch = as.data.frame(indexMM)
index_CAQ_out = as.data.frame(CAQ_index[-4,]) 
index_MM_out = indexMM

for (i in 1:nrow(index_CAQ_out)) {
  catch[,2] = index_CAQ_out[i,1]
  MM_CAQ_catch = rbind(MM_CAQ_catch,catch)
}

CAQ_out = CAQ[,-6]



MM_CAQ = as.data.frame(c())

for (i in 1:nrow(sym)) {
  print(i)
  catch3 = MM_CAQ_catch
  catch3[,3] = i
  catch = vector()
  
  for (n in 1:nrow(index_CAQ_out)) {
    catch = c(catch, t(sym[i,c(2:55)]))
  }
  
  
  catch2 = rep(t(CAQ_out[i,c(2:5)]), nrow(index_MM_out))
  
  
  catch3 = cbind(catch3, as.data.frame(catch), as.data.frame(catch2))
  
  MM_CAQ = rbind(MM_CAQ, catch3)
  
}




names(MM_CAQ) = c('MM','CHD','num','MM_type','CHD_type')

load('network.RData')

MM_CAQ_out = MM_CAQ



for (m in 1:nrow(MM_CAQ)) {
  print(m)
  MM_CAQ_out[m,6] = subset(symMM_CAQ_node, MM_CAQ_out$MM[m] == symMM_CAQ_node$lhs & MM_CAQ_out$CHD[m] == symMM_CAQ_node$rhs)$lift[1]
  MM_CAQ_out[m,7] = subset(symMM_CAQ_edge_list, MM_CAQ_out$MM[m] == symMM_CAQ_edge_list$source & MM_CAQ_out$CHD[m] == symMM_CAQ_edge_list$target)$correlation[1]
  
}

names(MM_CAQ_out) = c('MM','CHD','num','MM_type','CHD_type','lift','correlation')
MM_CAQ_out[is.na(MM_CAQ_out)] = 0

#save.image(file = "MM_CAQ_CNN.RData")


#library(openxlsx)
#write.xlsx(MM_CAQ_out, file = "MM_CAQ_out.xlsx", colNames = TRUE)
rm(list = ls())

load('MM_CAQ_CNN.RData')

num_type = as.data.frame(sym$num)
for (n in 1:nrow(num_type)) {
  print(n)
  catch = subset(MM_CAQ_out,MM_CAQ_out$num == n)
  num_type[n,2] = sum(catch$CHD_type)
}

CAQ_CNN = subset(num_type,num_type$V2 != 0)
NCAQ_CNN = subset(num_type, num_type$V2 == 0)

CAQ_CNN_data = MM_CAQ_out[which(MM_CAQ_out$num %in% CAQ_CNN$`sym$num`[1:425]),]
NCAQ_CNN_data = MM_CAQ_out[which(MM_CAQ_out$num %in% NCAQ_CNN$`sym$num`),]

library(openxlsx)
#write.xlsx(CAQ_CNN_data, file = "CAQ_CNN_data.xlsx", colNames = TRUE)
#write.xlsx(NCAQ_CNN_data, file = "NCAQ_CNN_data.xlsx", colNames = TRUE)

test_CNN = as.data.frame(c())

for (i in 1:(nrow(MM_CAQ_catch)*2)) {
  print(i)
  catch = MM_CAQ_catch
  catch$out = 0
  catch$out[i] = 1
  catch2 = MM_CAQ_catch
  catch2$out = 0
  catch = rbind(catch,catch2)
  test_CNN = rbind(test_CNN, catch)
  
}

library(openxlsx)
#write.xlsx(test_CNN, file = "test_CNN.xlsx", colNames = TRUE)




#save.image(file = "MM_CAQ_CNN$test.RData")



rm(list = ls())

load('MM_CAQ_CNN$test.RData')

CNN_out = read.csv('CHD_CAQ_meanPre.csv',header = T, stringsAsFactors = F)


catch$out = 1
catch2$out = 0
catch3 = rbind(catch,catch2)
names(catch3) = c('MM','CAQ','type')
CNN_out = CNN_out[,-1]
CNN_out = cbind(catch3, CNN_out)

names(MM_CAQ_catch) = c("lhs", "rhs")

rm(list=ls()[which(ls()!='CNN_out' & ls()!='MM_CAQ_catch')])

#save.image(file = "CNN_out.RData")






