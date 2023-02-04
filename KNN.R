rm(list = ls())

load('KNN_data.Rdata')


library(proxy)
library(dplyr)
library(kknn)
library(pROC)
library(modEvA)
library(class)

KNN_data$num = 1:nrow(KNN_data)
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

####KNN#####

KS = 6
#setInf = select(ZoutKNN,cut,weight,Euc,RWR,GoSimHY,GoSimHYSym)
setInf = ZoutKNN
library(kknn)
library(pROC)
library(modEvA)


library(class)

v = 10 #交叉次数

#分割原始数据（对原始数据均匀分成v重，以类别标签返回
#通过随机sample方法，如sample(5) --> 3 2 1 4 5,以sample产生的数据作为索引下标，提取cut类别标签数据。
#这样做到来随机抽取类别数据
grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
print(grps)
# 对每份数据分别运行ML函数 对随机数据的1/v作为测试集，剩下的作为训练集
# 对于lapply函数，后面要加上传递给function的参数data, cl, k. lapply的返回值为list
pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = fitted(kknn(cut ~ goSim + Euc + dij + CNN, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFKNN = setInfF
knn_roc = roc(setInfFKNN$cut, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfFKNN$cut,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
#aupr=AUC(obs=setInfFKNN$cut,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")



knn_roc$levels
setInfFKNN$knn = wh
setInfFKNN = setInfFKNN[order(setInfFKNN$num),]

ZoutKNN$knn = setInfFKNN$knn

#write.csv(setInfFKNN, 'KNNout.csv', quote = F, row.names = F, fileEncoding = 'GBK')


####KNNC1#####

KS = 6
#setInf = select(ZoutKNN,cut,weight,Euc,RWR,GoSimHY,GoSimHYSym)
setInf = ZoutKNN
library(kknn)
library(pROC)
library(modEvA)


library(class)

v = 10 #交叉次数

#分割原始数据（对原始数据均匀分成v重，以类别标签返回
#通过随机sample方法，如sample(5) --> 3 2 1 4 5,以sample产生的数据作为索引下标，提取cut类别标签数据。
#这样做到来随机抽取类别数据
grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
print(grps)
# 对每份数据分别运行ML函数 对随机数据的1/v作为测试集，剩下的作为训练集
# 对于lapply函数，后面要加上传递给function的参数data, cl, k. lapply的返回值为list
pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = fitted(kknn(cut ~ dij + CNN + Euc, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFKNN = setInfF
knn_roc = roc(setInfFKNN$cut, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfFKNN$cut,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
#aupr=AUC(obs=setInfFKNN$cut,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")



knn_roc$levels

knn_roc$levels
setInfFKNN$knnC1 = wh
setInfFKNN = setInfFKNN[order(setInfFKNN$num),]

ZoutKNN$knnC1 = setInfFKNN$knnC1

#write.csv(setInfFKNN, 'KNNout.csv', quote = F, row.names = F, fileEncoding = 'GBK')




####KNNC2#####

KS = 6
#setInf = select(ZoutKNN,cut,weight,Euc,RWR,GoSimHY,GoSimHYSym)
setInf = ZoutKNN
library(kknn)
library(pROC)
library(modEvA)


library(class)

v = 10 #交叉次数

#分割原始数据（对原始数据均匀分成v重，以类别标签返回
#通过随机sample方法，如sample(5) --> 3 2 1 4 5,以sample产生的数据作为索引下标，提取cut类别标签数据。
#这样做到来随机抽取类别数据
grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
print(grps)
# 对每份数据分别运行ML函数 对随机数据的1/v作为测试集，剩下的作为训练集
# 对于lapply函数，后面要加上传递给function的参数data, cl, k. lapply的返回值为list
pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = fitted(kknn(cut ~ goSim + num, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFKNN = setInfF
knn_roc = roc(setInfFKNN$cut, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfFKNN$cut,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
#aupr=AUC(obs=setInfFKNN$cut,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")



knn_roc$levels

knn_roc$levels
setInfFKNN$knnC2 = wh
setInfFKNN = setInfFKNN[order(setInfFKNN$num),]

ZoutKNN$knnC2 = setInfFKNN$knnC2

#write.csv(setInfFKNN, 'KNNout.csv', quote = F, row.names = F, fileEncoding = 'GBK')



######ROC#######


library(plotROC)
library(tidyverse)
library(ggplot2)
library(ggsci)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal
library("scales")
show_col(mypal)

aa = ZoutKNN


names(aa)

#####单条#####
basicplot <- ggplot(aa, aes(d = cut, m = knn )) + geom_roc(n.cuts = 0, labelsize = 5, labelround = 2)
basicplot
styledplot <- basicplot + style_roc()
styledplot
direct_label(basicplot, labels = "Biomarker", nudge_y = -.1) + style_roc()



########
aa$cut = -aa$cut
longtest <- melt_roc(aa, "cut", c("knn","knnC2","knnC1"))
longtest$name[which(longtest$name == "knnC1")] = "KNN Comparison1"
longtest$name[which(longtest$name == "knnC2")] = "KNN Comparison2"

longtest$name[which(longtest$name == "knn")] = "KNN"
head(longtest)
p2 = ggplot(longtest, aes(d = D, m = -M, color = name)) + 
  geom_roc(n.cuts = 0, labelsize = 5, labelround = 2) + 
  style_roc() +
  ggsci::scale_color_npg()
p2
auc<-calc_auc(p2)
head(auc)
p2 = p2+annotate("text",x = .70, y = .35, ## 注释text的位置
                 label = paste("AUC of KNN =", round(calc_auc(p2)$AUC[1], 3))) +
  annotate("text",x = .70, y = .25, ## 注释text的位置)
           label=paste("AUC of KNN Comparison1 =", round(calc_auc(p2)$AUC[2], 3)))+
  annotate("text",x = .70, y = .15, ## 注释text的位置)
           label=paste("AUC of KNN Comparison2 =", round(calc_auc(p2)$AUC[3], 3)))+
  ggtitle('ROC curve') +theme(plot.title = element_text(size = 17,hjust = 0.5))+
  theme(legend.title=element_blank())

p2


