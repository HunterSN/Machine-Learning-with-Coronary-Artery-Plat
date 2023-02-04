rm(list = ls())

library(dplyr)


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

MM_index = Gene_MM_CHD$MM[!duplicated(Gene_MM_CHD$MM)] %>% as.data.frame()
CHD_index = Gene_MM_CHD$CHD[!duplicated(Gene_MM_CHD$CHD)] %>% as.data.frame()

rm(list=ls()[which(ls()!='MM_index' & ls()!='CHD_index'& ls()!='symMM_CAQ')])

pinshu = as.data.frame(table(symMM_CAQ$index))




######MM饼图######

#MM
MM_pinshu = pinshu[which(pinshu$Var1 %in% MM_index[,1]),]

library(ggplot2)
library(ggsci)
library(plotROC)
library(tidyverse)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal
library("scales")
show_col(mypal)


gsfc = MM_pinshu
names(gsfc) = c('Symptom','pinshu')

gsfc = gsfc[order(gsfc$pinshu,decreasing = T),]
catch =  gsfc[which(gsfc$pinshu <= 343),]
catch2 = gsfc[which(gsfc$pinshu > 343),] 
catch3 = as.data.frame(c())
catch3[1,1] = 'Others'
catch3[1,2] = sum(catch$pinshu)
names(catch3) = c('Symptom','pinshu')
gsfc = rbind(catch2, catch3)


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )



myLabel = as.vector(gsfc$Symptom)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(gsfc$pinshu/ sum(gsfc$pinshu) * 100, 2), "%)", sep = "")   ## 用 round() 对结果保留两位小数



names(gsfc) = c('Symptom', 'Freq')

gsfc$Symptom = myLabel


gsfc$Symptom = factor(gsfc$Symptom, c("Chest Pain(28.55%)","Dyspnea(13.67%)","Xerostomia(6.61%)","Fatigue(6.08%)","Arteriosclerosis(5.59%)","Cardiac Arrhythmia(5.45%)", "Dizziness(5.11%)","Insomnia(5.06%)","Backache(3.22%)","Others(20.65%)"))


ggplot(data=gsfc, mapping=aes(x="Symptom",y=Freq,fill=Symptom))+
  geom_bar(stat="identity",width=0.5,position='stack',size=5)+
  coord_polar("y", start=0)+
  blank_theme+
  scale_fill_manual(values=c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2","#bebebe"))

#save(MM_pinshu, file = "MM_pinshu.RData")

######CAQ条形图######

#Coronary artery plaque
CAQ_pinshu = pinshu[which(pinshu$Var1 %in% CHD_index[,1]),]

library(ggplot2)
library(ggsci)
library(plotROC)
library(tidyverse)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal
library("scales")
show_col(mypal)


gsfc = CAQ_pinshu %>% as.matrix()


gsfc[1,1] = 'No coronary artery plaque'
gsfc = as.data.frame(gsfc)
gsfc$Freq = as.numeric(gsfc$Freq)
gsfc[1,2] = 1498 - gsfc[1,2] 
names(gsfc) = c('coronary artery plaque','Count')





gsfc$`coronary artery plaque` = factor(gsfc$`coronary artery plaque`, c("Calcified plaque","Non-calcified plaque","Mixed plaque","No coronary artery plaque"))

names(gsfc) = c("Class","Count")

ggplot(gsfc, aes(Class, Count))+ 
  geom_bar(stat="identity",aes(fill=Class),width=0.7,position='stack',size=5) +
  scale_fill_manual(values = c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2"))+
  theme(axis.title.x = element_blank(),
        text = element_text(size = 13),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))


#######网络图#######
rm(list = ls())

library(dplyr)


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

MM_index = Gene_MM_CHD$MM[!duplicated(Gene_MM_CHD$MM)] %>% as.data.frame()
CHD_index = Gene_MM_CHD$CHD[!duplicated(Gene_MM_CHD$CHD)] %>% as.data.frame()

rm(list=ls()[which(ls()!='MM_index' & ls()!='CHD_index'& ls()!='symMM_CAQ')])

pinshu = as.data.frame(table(symMM_CAQ$index))








load('network.RData')

symMM_CAQ_node = symMM_CAQ_node[,-2]
names(symMM_CAQ_node)[1:2] = c('source','target')
edge_net = symMM_CAQ_node[which((symMM_CAQ_node$source %in% MM_index[,1] &symMM_CAQ_node$target %in% CHD_index[,1])|(symMM_CAQ_node$target %in% MM_index[,1] &symMM_CAQ_node$source %in% CHD_index[,1])),]

#edge_net = symMM_CAQ_edge_list[which((symMM_CAQ_edge_list$source %in% MM_index[,1] &symMM_CAQ_edge_list$target %in% CHD_index[,1])|(symMM_CAQ_edge_list$target %in% MM_index[,1] &symMM_CAQ_edge_list$source %in% CHD_index[,1])),]


edge_net$source[which(edge_net$source == 'Coronary artery plaque')] = 'Non-calcified plaque'
edge_net$target[which(edge_net$target == 'Coronary artery plaque')] = 'Non-calcified plaque'

CHD_index = as.data.frame(CHD_index[2:4,])
CHD_index[,2] = 'Coronary artery plaque'
names(CHD_index) = c('node','type')
MM_index[,2] = 'Symptom'
names(MM_index) = c('node', 'type')



node_net = rbind(CHD_index,MM_index)


names(symMM_CAQ_node_list)[1] = "node"

node_net = merge(node_net,symMM_CAQ_node_list,all.x = T)


rm(list=ls()[which(ls()!='node_net' & ls()!='edge_net')])

KNN = read.csv('KNNout.csv',header = T)
KNN = KNN[,c(1,2,13)]

edge_net$KNN = 0

for (i in 1:nrow(edge_net)) {
  print(i)
  catch = subset(KNN,(KNN$lhs %in% edge_net$source[i] & KNN$rhs %in% edge_net$target[i])|(KNN$rhs %in% edge_net$source[i] & KNN$lhs %in% edge_net$target[i]))
  edge_net$KNN[i] = catch[1,3]
}




library(DelayedArray)
library(tidyverse)
library(igraph)


edges = edge_net

edges$KNN_type = 'Less than the cut-off value'
edges$KNN_type[which(edges$KNN >= 1)] = 'Greater the cut-off value'

vertices = node_net
vertices = vertices[which(vertices$degree != 'NA'),] 




names(edges)[1:2] = c('node1', 'node2')
names(vertices)[1:2] = c('node','type')




g <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE )


V(g)$color <- if_else(V(g)$type == "Symptom", "#E64B35B2", "#4DBBD5B2")
 
V(g)$size <- V(g)$degree*0.3

E(g)$weight <- abs(rnorm(n = ecount(g)))
E(g)$arrow.size <- .2
E(g)$edge.color <- "gray80"
graph_attr(g, "layout") <- layout_with_lgl



### 4.1 准备网络图布局数据
#?layout_in_circle # 帮助信息中，有其它布局函数。
layout1 <- layout_in_circle(g) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(g) # fr布局。
layout3 <- layout_on_grid(g) # grid布局。
head(layout1)

### 4.2 设置绘图颜色
#?rgb() # 可以将RGB值转换为16进制颜色名称。
#### 设置节点与分组背景色
color <- c(rgb(65,179,194,maxColorValue = 255),
           rgb(255,255,0,maxColorValue = 255),
           rgb(201,216,197,maxColorValue = 255))

names(color) <- unique(V(g)$type) # 将颜色以节点分类属性命名
V(g)$point.col <- color[match(V(g)$type,names(color))] # 设置节点颜色。
#names(color2) <- unique(V(g)$type) # 如果想要节点颜色与背景颜色不一致，则可以为节点单独设置一个颜色集。
#V(g)$point.col <- color2[match(V(g)$type,names(color2))]

#### 边颜色按照相关性正负设置
#E(g)$color <- ifelse(E(g)$linetype == "positive",rgb(255,215,0,maxColorValue = 255),"gray50")
E(g)$color <- ifelse(E(g)$KNN_type == 'Greater the cut-off value',"#F39B7FB2","#bebebe")

#pdf("Figure1_D.pdf",family = "Times",width = 10,height = 12)
#par(mar=c(5,2,1,2))


plot.igraph(g, layout=layout2,#更多参数设置信息?plot.igraph查看。
            ##节点颜色设置参数##
            vertex.color=V(g)$point.col,
            vertex.frame.color ="black",
            vertex.border=V(g)$point.col,
            ##节点大小设置参数##
            vertex.size=V(g)$n,
            ##节点标签设置参数##
            vertex.label=g$name,
            vertex.label.cex=1,
            #vertex.label.dist=0, # 标签距离节点中心的距离，0表示标签在节点中心。
            #vertex.label.degree = 0, # 标签对于节点的位置，0-右，pi-左，-pi/2-上，pi/2-下。
            vertex.label.col="black",
            ##链接属性参数以edge*开头##
            edge.arrow.size=0.5,
            edge.width=abs(E(g)$lift)*2,
            edge.curved = TRUE
)
##设置图例，与plot.igraph()函数一起运行##
legend(
  title = "Type",
  list(x = min(layout1[,1])-0.2,
       y = min(layout1[,2])-0.17), # 图例的位置需要根据自己的数据进行调整，后面需要使用AI手动调整。
  legend = c(unique(V(g)$type)),
  fill = color,
  #pch=1
)

legend(
  title = "Lift",
  list(x = min(layout1[,1])+0.4,
       y = min(layout1[,2])-0.17),
  legend = c(1,2,3,4),
  col = "black",
  lty=1,
  lwd=c(1,2,3,4)*2,
)


legend(
  title = 'KNN type',
  list(x = min(layout1[,1])+0.6,
       y = min(layout1[,2])-0.17),
  legend = c("Greater the cut-off value","Less than the cut-off value"),
  col = c("#F39B7FB2","#bebebe"),
  lty=1,
  lwd=4
)
legend(
  title = "Degree",
  list(x = min(layout1[,1])+1.4,
       y = min(layout1[,2])-0.17),
  legend = c(1,seq(0,8,2)[-1]),# max(V(g)$degree)
  col = "black",
  pch=1,
  pt.lwd=1,
  pt.cex=c(1,seq(0,8,2)[-1])
)

#dev.off()





#write.csv(edges, 'edges.csv', quote = F, row.names = F, fileEncoding = 'GBK')




