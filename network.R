rm(list = ls())

library(Matrix)
library(arules) #apriori
library(arulesViz)
library(igraph)
library(car)

load('sym&CAQ_Process.RData')
rm(list=ls()[which(ls()!='symMM' & ls() != 'symMM_CAQ')])

######获取节点信息#####
#MM
a = symMM#读取数据
b_apriori = split(a$index,a$num)#按处方分类组成list
c_apriori = apriori(as(b_apriori,'transactions'),parameter = list(supp = 0,conf = 0,maxlen = 2,minlen = 2)) #apriori计算
summary(c_apriori)
output =inspect(c_apriori)

shai = output
shai$lhs =  gsub('[{]',"",shai$lhs)
shai$lhs =  gsub('[}]',"",shai$lhs)

shai$rhs =  gsub('[{]',"",shai$rhs)
shai$rhs =  gsub('[}]',"",shai$rhs)

symMM_node = shai[which(shai$count != 0),]

#MM&CAQ

a = symMM_CAQ#读取数据
b_apriori = split(a$index,a$num)#按处方分类组成list
c_apriori = apriori(as(b_apriori,'transactions'),parameter = list(supp = 0,conf = 0,maxlen = 2,minlen = 2)) #apriori计算
summary(c_apriori)
output =inspect(c_apriori)

shai = output
shai$lhs =  gsub('[{]',"",shai$lhs)
shai$lhs =  gsub('[}]',"",shai$lhs)

shai$rhs =  gsub('[{]',"",shai$rhs)
shai$rhs =  gsub('[}]',"",shai$rhs)

symMM_CAQ_node = shai[which(shai$count != 0),]
rm(list=ls()[which(ls()!='symMM' & ls() != 'symMM_CAQ' & ls() != 'symMM_node' & ls() != 'symMM_CAQ_node')])



#######构造关系网络(简单版)#######

data = symMM_CAQ_node[,c(1,3)]

dir=F

rem.multi=T

labels<-union(unique(data[,1]),unique(data[,2]))
ids<-1:length(labels);names(ids)<-labels
from<-as.character(data[,1]);to<-as.character(data[,2])
edges<-matrix(c(ids[from],ids[to]),nc=2)
g<-graph.empty(directed = dir)
g<-add.vertices(g,length(labels))
V(g)$label=labels
g<-add.edges(g,t(edges))
if (rem.multi){
  E(g)$weight<-count.multiple(g)
  g<-simplify(g,remove.multiple = TRUE,
              remove.loops = TRUE,edge.attr.comb = "mean")
}
g





init.igraph<-function(data,dir=F,rem.multi=T){
  labels<-union(unique(data[,1]),unique(data[,2]))
  ids<-1:length(labels);names(ids)<-labels
  from<-as.character(data[,1]);to<-as.character(data[,2])
  edges<-matrix(c(ids[from],ids[to]),nc=2)
  g<-graph.empty(directed = dir)
  g<-add.vertices(g,length(labels))
  V(g)$label=labels
  g<-add.edges(g,t(edges))
  if (rem.multi){
    E(g)$weight<-count.multiple(g)
    g<-simplify(g,remove.multiple = TRUE,
                remove.loops = TRUE,edge.attr.comb = "mean")
  }
  g
}

a =  init.igraph(data)


#####建立矩阵#######
#MM
a = symMM_node
index = rbind(as.matrix(a$lhs), as.matrix(a$rhs))
index = as.matrix(index[!duplicated(index)])
#index = as.matrix(as.data.frame(table(GeneDiseaseRelations_PPI$mim_morbid))[,1])


juzhen = matrix(1:nrow(index)^2, nrow = nrow(index), dimnames = list(index, index))
juzhen = as.data.frame(juzhen)

for (lie in 1:ncol(juzhen)) {
  print(lie)
  for (hang in 1:nrow(juzhen)) {
    mm = subset(a, (a$lhs == names(juzhen)[lie] & a$rhs == colnames(juzhen)[hang]) | (a$rhs == names(juzhen)[lie] & a$lhs == colnames(juzhen)[hang]))
    juzhen[hang, lie] = mm$count[1]
  }
}

juzhen[is.na(juzhen)] = 0

symMM_matrix = juzhen

#MM&CAQ


a = symMM_CAQ_node
index = rbind(as.matrix(a$lhs), as.matrix(a$rhs))
index = as.matrix(index[!duplicated(index)])

juzhen = matrix(1:nrow(index)^2, nrow = nrow(index), dimnames = list(index, index))
juzhen = as.data.frame(juzhen)

for (lie in 1:ncol(juzhen)) {
  print(lie)
  for (hang in 1:nrow(juzhen)) {
    mm = subset(a, (a$lhs == names(juzhen)[lie] & a$rhs == colnames(juzhen)[hang]) | (a$rhs == names(juzhen)[lie] & a$lhs == colnames(juzhen)[hang]))
    juzhen[hang, lie] = mm$count[1]
  }
}

juzhen[is.na(juzhen)] = 0

symMM_CAQ_matrix = juzhen

rm(list =  c('a','index','juzhen','mm','hang','lie'))

#rm(list=ls()[which(ls()!='symMM_CAQ_matrix' & ls()!='symMM_matrix')])
#save.image(file = "MM&CAQ_matrix.RData")

#######网络分析############
library(igraph)

#输入数据示例，邻接矩阵
#这是一个微生物相关性网络，数值表示了微生物 OTU 之间的相关系数（正负相关以及相关性的大小）
adjacency_weight <- symMM_CAQ_matrix
#adjacency_weight <-symMM_matrix
head(adjacency_weight)[1:6]    #邻接矩阵类型的网络文件

#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络
igraph = graph_from_adjacency_matrix(as.matrix(adjacency_weight), mode = 'undirected', weighted = TRUE, diag = FALSE)
igraph    #igraph 的邻接列表

#上述将相关系数转化为边的权重
#由于相关系数有负值，而权重正常都应为正值，所以对由相关系数得到的权重取绝对值
#新生成一列边属性记录相关系数
E(igraph)$corr <- E(igraph)$weight
E(igraph)$weight <- abs(E(igraph)$weight)

##节点特征
#节点数量
length(V(igraph)$name)
#或
vcount(igraph)

#节点度（Degree）
#由于本示例是个无向网络，故无出度和入度之分
V(igraph)$degree <- degree(igraph)
V(igraph)$degree

#查看度分布
#可观察到微生物相关网络通常服从幂律分布，这个下节再讲怎样通过计算验证
degree_dist <- degree.distribution(igraph)[-1]
degree_num <- 1:max(V(igraph)$degree)

par(mfrow = c(1, 2))
hist(V(igraph)$degree, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_num, degree_dist, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-intensity', main = 'Log-log degree distribution')

#查看节点度与其“邻居”的平均度的关系
#微生物网络中高度值的节点更倾向连接在一起，是普遍现象吗？
neighbor_degree <- graph.knn(igraph, V(igraph))$knn
plot(V(igraph)$degree, neighbor_degree, log = 'xy', 
     xlab = 'Log degree', ylab = 'Log average neighbor degree')

#加权度（Weighted degree）
V(igraph)$weight_degree <- strength(igraph)
V(igraph)$weight_degree

#接近中心性（Closeness centrality）
V(igraph)$closeness_centrality <- closeness(igraph)
V(igraph)$closeness_centrality

#介数中心性（Betweenness centrality）
V(igraph)$betweenness_centrality <- betweenness(igraph)
V(igraph)$betweenness_centrality

#特征向量中心性（Eigenvector centrality）
V(igraph)$eigenvector_centrality <- evcent(igraph)$vector
V(igraph)$eigenvector_centrality

#探索三种描述节点中心性的特征的关系
library(car)

scatter3d(V(igraph)$closeness_centrality, V(igraph)$betweenness_centrality, V(igraph)$eigenvector_centrality, 
          xlab =  'Closeness centrality', ylab = 'Betweenness centrality', zlab = 'Eigenvector centrality', 
          surface = FALSE)

#探索节点度和节点中心性的关系，如与特征向量中心性的关系
plot(V(igraph)$degree, V(igraph)$eigenvector_centrality, 
     xlab = 'Degree', ylab = 'Eigenvector centrality')

#输出列表
node_list <- data.frame(
  node_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  weight_degree = V(igraph)$weight_degree, 
  closeness_centrality = V(igraph)$closeness_centrality, 
  betweenness_centrality = V(igraph)$betweenness_centrality, 
  eigenvector_centrality = V(igraph)$eigenvector_centrality)

head(node_list)
#write.table(node_list, 'node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##边特征
#边的数量
ecount(igraph)

#权重（Weighted），已在数据读入时转化获得
E(igraph)$weight

#边介数中心性（Edge betweenness centrality）
E(igraph)$betweenness_centrality <- edge.betweenness(igraph)
E(igraph)$betweenness_centrality

#输出列表
edge <- data.frame(as_edgelist(igraph))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(igraph)$weight,
  correlation = E(igraph)$corr, 
  betweenness_centrality = E(igraph)$betweenness_centrality
)
head(edge_list)

#write.table(edge_list, 'edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)


symMM_CAQ_edge_list = edge_list
symMM_CAQ_node_list = node_list


rm(list =  c('adjacency_weight', 'edge', 'igraph', 'degree_dist', 'degree_num', 'neighbor_degree', 'edge_list', 'node_list'))


#save.image(file = "network.RData")









#save.image(file = "zgxfTapShaiApriori.RData")

#library(openxlsx)
#write.xlsx(shai, file = "zgxfZongApriori.xlsx", colNames = TRUE)

#load('knnData.RData')
#rm(list=ls()[which(ls()!='ZoutKNN')])

#tap筛节点
#taps = read.csv('taps.csv',header = F,stringsAsFactors = F)
#aa = subset(shai, (is.element(shai$lhs,taps[,1])) & (is.element(shai$rhs,taps[,1])))


#catch[is.na(catch)] = 0
