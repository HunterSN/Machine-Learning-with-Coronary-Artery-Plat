rm(list = ls())

load('sym&CAQ_Process.RData')
rm(list=ls()[which(ls()!='symMM_CAQ')])

trait = read.csv('source2.csv',header = T)


load('MM_CAQ_goSim.Rdata')
genes = rbind(as.matrix(Gene_CAQ_group$labels),as.matrix(Gene_Ca_group$labels),as.matrix(Gene_NCa_group$labels),as.matrix(Gene_Mix_group$labels))
genes = as.data.frame(genes[!duplicated(genes)])
names(genes) = c('substanceBXH')

Gene_CAQ_group$CHD = 'Coronary artery plaque'
Gene_Ca_group$CHD = 'Calcified plaque'
Gene_NCa_group$CHD = 'Non-calcified plaque'
Gene_Mix_group$CHD = 'Mixed plaque'


Gene_MM_CHD =  rbind(as.matrix(Gene_CAQ_group),as.matrix(Gene_Ca_group),as.matrix(Gene_NCa_group),as.matrix(Gene_Mix_group)) %>% as.data.frame()

Z_index = rbind(as.matrix(Gene_MM_CHD$MM), as.matrix(Gene_MM_CHD$CHD)) 
Z_index = Z_index[!duplicated(Z_index)] %>% as.data.frame()

symMM_CAQ = symMM_CAQ[which(symMM_CAQ$index %in% Z_index[,1]),]

P_count = as.data.frame(symMM_CAQ$num[!duplicated(symMM_CAQ$num)])

MM_index = Gene_MM_CHD$MM[!duplicated(Gene_MM_CHD$MM)]
CHD_index = Gene_MM_CHD$CHD[!duplicated(Gene_MM_CHD$CHD)]





for (i in 1:nrow(P_count)) {
  print(i)
  catch = as.data.frame(c())
  catch2 = symMM_CAQ[which(symMM_CAQ$num == P_count[i,1]),]
  catch3 = subset(catch2,catch2$index %in% MM_index)
  catch4 = subset(catch2,catch2$index %in% CHD_index)
  
  catch5 = subset(Gene_MM_CHD, (Gene_MM_CHD$MM %in% catch3$index) & (Gene_MM_CHD$CHD %in% catch4$index))
  
  
  for (m in 1:nrow(genes)) {
    catch[m,1] = subset(catch5,catch5$labels == genes[m,1]) %>% nrow()
  }
  names(catch) = c(paste('P',i, sep = '_'))
  genes = cbind(genes,catch)
  
}





trait = trait[which(trait$num %in% P_count[,1]),]
for (i in 1:nrow(trait)) {
  trait[i,1] = paste('P',trait[i,1],sep = '_')
}


library(mice)
library(DMwR2)

D_NA = function(x){
  x = dplyr::mutate_all(x,as.numeric)
  miceMod <- mice(x[, !names(x) %in% c("CASEID")], method="rf")  # 基于随机森林模型进行mice插值
  miceOutput <- complete(miceMod)  # 生成完整数据
  return(miceOutput)
}


trait2 = D_NA(trait)


trait2$num = trait$num

#rm(list=ls()[which(ls()!='trait2' & ls()!='genes')])
#save.image(file = "WGCNA_data.RData")



