rm(list = ls())

load('PPIN_SYMBOL.RData')

PPIN_SYMBOL = subset(PPIN_SYMBOL, PPIN_SYMBOL$combined_score >= 400)

PPIN = PPIN_SYMBOL[,c(15:16)]

rm(list = c('PPIN_SYMBOL', 'proInfo'))

MM_target = read.csv('MM_target.csv',header = T, stringsAsFactors = F)
gene = as.matrix(MM_target$Gene.symbol)
gene = as.matrix(gene[!duplicated(gene)])

PPIN_Out = subset(PPIN, PPIN$source %in% gene & PPIN$target %in% gene)
PPIN_gene = rbind(as.matrix(PPIN_Out$source), as.matrix(PPIN_Out$target))
PPIN_gene = as.matrix(PPIN_gene[!duplicated(PPIN_gene)])

MM_target_Net = cbind(as.matrix(MM_target$Gene.symbol), as.matrix(MM_target$MM))

load('network.RData')

MM_CAQ_Net = cbind(as.matrix(symMM_CAQ_node$lhs), as.matrix(symMM_CAQ_node$rhs))


MM_target_CAQ_Net = rbind(MM_target_Net, MM_CAQ_Net)

PPIN_Out = as.matrix(PPIN_Out)

PPIN_MM_target_CAQ_Net = rbind(PPIN_Out, MM_target_CAQ_Net)


#####建立矩阵#######



a = as.data.frame(PPIN_MM_target_CAQ_Net)

index = as.matrix(names(symMM_CAQ_matrix))

#tap = read.csv('CAQ_index.csv',header = T)
#tap = as.matrix(tap)



juzhen = matrix(1:nrow(index)^2, nrow = nrow(index), dimnames = list(index, index))
juzhen = as.data.frame(juzhen)

for (lie in 1:ncol(juzhen)) {
  print(lie)
  for (hang in 1:nrow(juzhen)) {
    mm = subset(a, (a$source == names(juzhen)[lie] & a$target == colnames(juzhen)[hang]) | (a$target == names(juzhen)[lie] & a$source == colnames(juzhen)[hang]))
    juzhen[hang, lie] = nrow(mm)
  }
}

juzhen[is.na(juzhen)] = 0

symMM_CAQ_PPIN_matrix = juzhen

rm(list =  c('a','index','juzhen','mm','hang','lie'))

symMM_CAQ_PPIN_matrix[symMM_CAQ_PPIN_matrix == 2] = 1

#######网络分析############
library(igraph)

#输入数据示例，邻接矩阵
#这是一个微生物相关性网络，数值表示了微生物 OTU 之间的相关系数（正负相关以及相关性的大小）
adjacency_weight <- symMM_CAQ_PPIN_matrix
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


symMM_CAQ_PPIN_edge_list = edge_list
symMM_CAQ_PPIN_node_list = node_list


rm(list =  c('adjacency_weight', 'edge', 'igraph', 'degree_dist', 'degree_num', 'neighbor_degree', 'edge_list', 'node_list'))


#save.image(file = "PPIN_network.RData")



load('PPIN_network.RData')








########STRING_online########

library(tidyverse)
library(clusterProfiler) # Y叔的包有没有，这个其实只是为了ID转换
library(org.Hs.eg.db)  #小鼠的话，把Hs改成Mm

library(STRINGdb)
library(igraph)
library(ggraph)
# 创建STRINGdb对象
string_db <- STRINGdb$new( version="11.5", species=9606, 
                           score_threshold=400, input_directory="")





# Y叔的clusterProfiler将Gene Symbol转换为Entrez ID
gene <- gene %>% bitr(fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db", 
                      drop = T)



#使用map函数用于将基因匹配到STRING数据库的ID，map函数的帮助信息可以查看STRINGdb$help("map")。
#然后plot_network绘图即可。中途需要在线下载数据，慢慢等吧。


data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                      removeUnmappedRows = TRUE)


string_db$plot_network(data_mapped$STRING_id )















