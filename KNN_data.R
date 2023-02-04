rm(list = ls())

load('MM&CAQ_matrix.RData')
load('CNN_out.RData')



########距离######
nutrient = symMM_CAQ_matrix
herbCoup = MM_CAQ_catch
####欧式距离####
dEuc<-dist(nutrient, method = 'euclidean')
dEuc = as.matrix(dEuc)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Euc[i] = as.character(subset(dEuc, rownames(dEuc) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}

#####马氏距离######
dMan<-dist(nutrient, method = 'manhattan')
dMan = as.matrix(dMan)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Man[i] = as.character(subset(dMan, rownames(dMan) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}


#####堪培拉距离######
dCan<-dist(nutrient, method = 'canberra')
dCan = as.matrix(dCan)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Can[i] = as.character(subset(dCan, rownames(dCan) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}

#####二元变量的距离######
dBin<-dist(nutrient, method = 'binary')
dBin = as.matrix(dBin)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Bin[i] = as.character(subset(dBin, rownames(dBin) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}

#write.csv(herbCoup,"HNherbCoupOut.csv",quote = F,row.names = F,fileEncoding = 'GBK')
MM_CAQ_dis = herbCoup
rm(list=ls()[which(ls()!='MM_CAQ_dis')])
#save.image(file = "MM_CAQ_dis.RData")


#####GoSim######

rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)


load('CNN_out.RData')
load('MM_CAQ_dis.RData')

load('RWR_PPIN_MM_target_CAQ_seedCAQ.RData')
load('RWR_PPIN_MM_target_CAQ_seed_Calcified_plaque.RData')
load('RWR_PPIN_MM_target_CAQ_seed_Non_calcified_plaque.RData')
load('RWR_PPIN_MM_target_CAQ_seed_Mixed_plaque.RData')

symbol_index = read.csv('symbol_Index.csv',header = T, stringsAsFactors = F)
MM_index = as.data.frame(symbol_index$MM[which(symbol_index$MM != '无')])
CAQ_index = read.csv('CAQ_index.csv',header = T)

MM_target = read.csv('MM_target.csv',header = T, stringsAsFactors = F)
MM_target = MM_target[,c(2,13)]

#gene

Gene_CAQ_group = RWR_PPIN_MM_target_CAQ_seedCAQ[!(RWR_PPIN_MM_target_CAQ_seedCAQ$labels %in% CAQ_index[,1]),]
Gene_CAQ_group = Gene_CAQ_group[!(Gene_CAQ_group$labels %in% MM_index[,1]),]



Gene_Ca_group = RWR_PPIN_MM_target_CAQ_seed_Calcified_plaque[!(RWR_PPIN_MM_target_CAQ_seed_Calcified_plaque$labels %in% CAQ_index[,1]),]
Gene_Ca_group = Gene_Ca_group[!(Gene_Ca_group$labels %in% MM_index[,1]),]



Gene_NCa_group = RWR_PPIN_MM_target_CAQ_seed_Non_calcified_plaque[!(RWR_PPIN_MM_target_CAQ_seed_Non_calcified_plaque$labels %in% CAQ_index[,1]),]
Gene_NCa_group = Gene_NCa_group[!(Gene_NCa_group$labels %in% MM_index[,1]),]

Gene_Mix_group = RWR_PPIN_MM_target_CAQ_seed_Mixed_plaque[!(RWR_PPIN_MM_target_CAQ_seed_Mixed_plaque$labels %in% CAQ_index[,1]),]
Gene_Mix_group = Gene_Mix_group[!(Gene_Mix_group$labels %in% MM_index[,1]),]


#Gene_CAQ_group = Gene_CAQ_group[order(Gene_CAQ_group$Score, decreasing = T),]


#CAQ
Gene_CAQ_group = Gene_CAQ_group[which(Gene_CAQ_group$Score >= quantile(Gene_CAQ_group$Score)[4]),]
Gene_CAQ_group = Gene_CAQ_group[which(Gene_CAQ_group$Score >= quantile(Gene_CAQ_group$Score)[4]),]
Gene_CAQ_group = Gene_CAQ_group[which(Gene_CAQ_group$Score >= quantile(Gene_CAQ_group$Score)[4]),]
Gene_CAQ_group = Gene_CAQ_group[which(Gene_CAQ_group$Score >= quantile(Gene_CAQ_group$Score)[4]),]

catch = MM_target
catch = MM_target[which(MM_target$Gene.symbol %in% Gene_CAQ_group$labels),]
names(catch) = c('labels','MM')

Gene_CAQ_group = merge(Gene_CAQ_group, catch)
table(Gene_CAQ_group$MM)

#Ca
Gene_Ca_group = Gene_Ca_group[which(Gene_Ca_group$Score >= quantile(Gene_Ca_group$Score)[4]),]
Gene_Ca_group = Gene_Ca_group[which(Gene_Ca_group$Score >= quantile(Gene_Ca_group$Score)[4]),]
Gene_Ca_group = Gene_Ca_group[which(Gene_Ca_group$Score >= quantile(Gene_Ca_group$Score)[4]),]
Gene_Ca_group = Gene_Ca_group[which(Gene_Ca_group$Score >= quantile(Gene_Ca_group$Score)[4]),]

catch = MM_target
catch = MM_target[which(MM_target$Gene.symbol %in% Gene_Ca_group$labels),]
names(catch) = c('labels','MM')

Gene_Ca_group = merge(Gene_Ca_group, catch)
table(Gene_Ca_group$MM)



#NCa
Gene_NCa_group = Gene_NCa_group[which(Gene_NCa_group$Score >= quantile(Gene_NCa_group$Score)[4]),]
Gene_NCa_group = Gene_NCa_group[which(Gene_NCa_group$Score >= quantile(Gene_NCa_group$Score)[4]),]
Gene_NCa_group = Gene_NCa_group[which(Gene_NCa_group$Score >= quantile(Gene_NCa_group$Score)[4]),]
Gene_NCa_group = Gene_NCa_group[which(Gene_NCa_group$Score >= quantile(Gene_NCa_group$Score)[4]),]

catch = MM_target
catch = MM_target[which(MM_target$Gene.symbol %in% Gene_NCa_group$labels),]
names(catch) = c('labels','MM')

Gene_NCa_group = merge(Gene_NCa_group, catch)
table(Gene_NCa_group$MM)

#Mix
Gene_Mix_group = Gene_Mix_group[which(Gene_Mix_group$Score >= quantile(Gene_Mix_group$Score)[4]),]
Gene_Mix_group = Gene_Mix_group[which(Gene_Mix_group$Score >= quantile(Gene_Mix_group$Score)[4]),]
Gene_Mix_group = Gene_Mix_group[which(Gene_Mix_group$Score >= quantile(Gene_Mix_group$Score)[4]),]
Gene_Mix_group = Gene_Mix_group[which(Gene_Mix_group$Score >= quantile(Gene_Mix_group$Score)[4]),]

catch = MM_target
catch = MM_target[which(MM_target$Gene.symbol %in% Gene_Mix_group$labels),]
names(catch) = c('labels','MM')

Gene_Mix_group = merge(Gene_Mix_group, catch)
table(Gene_Mix_group$MM)
#####1#########
#Gene_MM = rbind(Gene_Ca_group,Gene_CAQ_group,Gene_Mix_group,Gene_NCa_group)
#Gene_MM = Gene_MM[!duplicated(Gene_MM),]
#Gene_MM = Gene_MM[,c(1,4)]
#Gene_MM = Gene_MM[!duplicated(Gene_MM),]
#rm(list=ls()[which(ls()!='Gene_MM')])
#save.image(file = 'Gene_MM.Rdata')


goSim = function(gene){
  keytypes(org.Hs.eg.db)
  gene.df <- bitr(gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                  toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                  OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
  head(gene.df)
  d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
  genesim = gene.df[,-2]
  genesim = genesim[!duplicated(genesim$SYMBOL),]
  num = unique(as.vector(genesim$ENTREZID))
  genetolist = t(combn(num,2)) 
  rn = nrow(genetolist)
  genetolist = as.data.frame(genetolist)
  for (l in (1:rn)) {
    catch = geneSim(genetolist[l,1],genetolist[l,2],semData=d, measure = "Wang", drop = "IEA", combine = "BMA")
    genetolist[l,3] = catch[1]
    #print(l)
  }
  names(genetolist) = c("geneA","geneB","goSim")
  genetolist
}

#CAQ
dataa = Gene_CAQ_group
CAQlist = as.matrix(dataa$MM)
CAQlist = as.matrix(CAQlist[!duplicated(CAQlist)])

gosim_CAQ_out = as.data.frame(CAQlist)
gosim_CAQ_out$goSim = 0
gosim_CAQ_out$CHD = 'Coronary artery plaque'

for (i in 1:nrow(CAQlist)) {
  print(i)
  genelist = dataa[which(dataa$MM == CAQlist[i,1]),]
  gene = genelist[,1]
  catch = goSim(gene)
  gosim_CAQ_out[i,2] = sum(catch$goSim)
  
}



#Ca

dataa = Gene_Ca_group
Calist = as.matrix(dataa$MM)
Calist = as.matrix(Calist[!duplicated(Calist)])

gosim_Ca_out = as.data.frame(Calist)
gosim_Ca_out$goSim = 0
gosim_Ca_out$CHD = 'Calcified plaque'

for (i in 1:nrow(Calist)) {
  print(i)
  genelist = dataa[which(dataa$MM == Calist[i,1]),]
  gene = genelist[,1]
  catch = goSim(gene)
  gosim_Ca_out[i,2] = sum(catch$goSim)
  
}



#NCa

dataa = Gene_NCa_group
NCalist = as.matrix(dataa$MM)
NCalist = as.matrix(NCalist[!duplicated(NCalist)])

gosim_NCa_out = as.data.frame(NCalist)
gosim_NCa_out$goSim = 0
gosim_NCa_out$CHD = 'Non-calcified plaque'

for (i in 1:nrow(NCalist)) {
  print(i)
  genelist = dataa[which(dataa$MM == NCalist[i,1]),]
  gene = genelist[,1]
  catch = goSim(gene)
  gosim_NCa_out[i,2] = sum(catch$goSim)
  
}


#Mix

dataa = Gene_Mix_group
Mixlist = as.matrix(dataa$MM)
Mixlist = as.matrix(Mixlist[!duplicated(Mixlist)])

gosim_Mix_out = as.data.frame(Mixlist)
gosim_Mix_out$goSim = 0
gosim_Mix_out$CHD = 'Mixed plaque'

for (i in 1:nrow(Mixlist)) {
  print(i)
  genelist = dataa[which(dataa$MM == Mixlist[i,1]),]
  gene = genelist[,1]
  catch = goSim(gene)
  gosim_Mix_out[i,2] = sum(catch$goSim)
  
}

goSim_out = rbind(gosim_CAQ_out, gosim_Ca_out,gosim_NCa_out,gosim_Mix_out)
names(goSim_out) = c("lhs","goSim", "rhs")

goSim_out = merge(MM_CAQ_catch ,goSim_out, all.x = T)
goSim_out[is.na(goSim_out)] = 0

#save.image(file = 'MM_CAQ_goSim.Rdata')
####RWR_process######



RWR_MM_CAQ_group = RWR_PPIN_MM_target_CAQ_seedCAQ[(RWR_PPIN_MM_target_CAQ_seedCAQ$labels %in% MM_index[,1]),]
RWR_MM_CAQ_group$CHD = 'Coronary artery plaque'

RWR_MM_Ca_group = RWR_PPIN_MM_target_CAQ_seed_Calcified_plaque[(RWR_PPIN_MM_target_CAQ_seed_Calcified_plaque$labels %in% MM_index[,1]),]
RWR_MM_Ca_group$CHD = 'Calcified plaque'

RWR_MM_NCa_group = RWR_PPIN_MM_target_CAQ_seed_Non_calcified_plaque[(RWR_PPIN_MM_target_CAQ_seed_Non_calcified_plaque$labels %in% MM_index[,1]),]
RWR_MM_NCa_group$CHD = 'Non-calcified plaque'

RWR_MM_Mix_group = RWR_PPIN_MM_target_CAQ_seed_Mixed_plaque[(RWR_PPIN_MM_target_CAQ_seed_Mixed_plaque$labels %in% MM_index[,1]),]
RWR_MM_Mix_group$CHD = 'Mixed plaque'


RWR_MM_CAQ = rbind(RWR_MM_CAQ_group, RWR_MM_Ca_group,RWR_MM_NCa_group,RWR_MM_Mix_group)
RWR_MM_CAQ = RWR_MM_CAQ[,-1]
names(RWR_MM_CAQ) = c("RWR","lhs", "rhs")

RWR_MM_CAQ_out = as.data.frame(MM_CAQ_catch)


RWR_MM_CAQ_out = merge(RWR_MM_CAQ_out, RWR_MM_CAQ, all.x = T)


RWR_MM_CAQ_out[is.na(RWR_MM_CAQ_out)] = 0



#rm(list=ls()[which(ls()!='goSim_out' & ls()!='RWR_MM_CAQ_out')])


#save.image(file = 'MM_CAQ_goSim&RWR.Rdata')

######dijkstra##########
rm(list = ls())

load('MM&CAQ_matrix.RData')
load('CNN_out.RData')

load('MM_CAQ_dijkstra.RData')


catch = MM_CAQ_catch
catch$dij = 0

for (i in 1:nrow(MM_CAQ_catch)) {
  print(i)
  catch$dij[i] = as.character(subset(symMM_CAQ_matrix, rownames(symMM_CAQ_matrix) == catch[i,1], select = as.character(catch[i,2])))
  
}

catch$dij[which(catch$dij == 'numeric(0)')] = 0


MM_CAQ_dij = catch

#rm(list=ls()[which(ls()!='MM_CAQ_dij')])
#save.image(file = 'MM_CAQ_dij.Rdata')


######合并#########

rm(list = ls())

load('MM_CAQ_dis.RData')
load('MM_CAQ_dij.Rdata')
load('MM_CAQ_goSim&RWR.Rdata')

load('CNN_out.RData')







a = as.data.frame(table(CNN_out$X0))

CNN_out$type[which(CNN_out$X0 == a$Var1[2])] = 0
CNN_out = CNN_out[1:216,1:3]
names(CNN_out) = c("lhs", "rhs", "CNN")

KNN_data = merge(CNN_out, goSim_out)
KNN_data = merge(KNN_data, MM_CAQ_dis)
KNN_data = merge(KNN_data, MM_CAQ_dij)
KNN_data = merge(KNN_data, RWR_MM_CAQ_out)

KNN_data = KNN_data[!duplicated(KNN_data),]

#rm(list=ls()[which(ls()!='KNN_data')])
#save.image(file = 'KNN_data.Rdata')



