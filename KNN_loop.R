rm(list = ls())

load('KNN_data.Rdata')


library(proxy)
library(dplyr)
library(kknn)
library(pROC)
library(modEvA)
library(class)


KNN_data$Euc = as.numeric(KNN_data$Euc)
KNN_data$Man = as.numeric(KNN_data$Man)
KNN_data$Can = as.numeric(KNN_data$Can)
KNN_data$Bin = as.numeric(KNN_data$Bin)
KNN_data$CNN = as.numeric(KNN_data$CNN)
KNN_data$RWR = as.numeric(KNN_data$RWR)
KNN_data$dij = as.numeric(KNN_data$dij)

ZoutKNN = KNN_data


#确定分割点
ZoutKNN$cut = 1
ZoutKNN$cut[which(ZoutKNN$RWR == 0)] = 0
rownames(ZoutKNN) = 1:nrow(ZoutKNN)


set.seed(1234)


par(mfrow = c(1,1))

####KNN循环#####
setInf = ZoutKNN
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(cut ~ goSim + Euc + dij + CNN, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  wh = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(wh))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$cut, as.numeric(wh),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$cut, as.numeric(wh))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)

