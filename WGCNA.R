rm(list = ls())

library(WGCNA)
#这个设置很重要，不要忽略，允许R语言以最大线程运行
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

######数据导入、清洗及预处理#######

load('WGCNA_data.RData')


# 在读入数据时，遇到字符串之后，不将其转换为factors，仍然保留为字符串格式
options(stringsAsFactors = FALSE)
# 导入示例数据，这里填写自己存放表达矩阵的路径
femData = genes	# femData代表该文件
# 查看数据
dim(femData)	# dim查看矩阵形状
names(femData)	# names查看列标，即每一列的标题



# 提取样本-基因表达矩阵
datExpr0 = as.data.frame(t(femData[, -c(1)]))
# 删除femData矩阵第1到8列，再转置，再变为数据格式，将所得用datExpr0表示
# 第一列是基因名，也删除了，即要以所有列标都是样本名的形式做转置
names(datExpr0) = femData$substanceBXH
# 转置后，行标变成了样本名，然后重新加入基因名作为列标
# femData$substanceBXH表示在femData矩阵里面依次取substanceBXH列的值
# names表示列标，即将这些值（基因名）作为datExpr0的列标
rownames(datExpr0) = names(femData)[-c(1)]
# 将femData列标的第1到8个删除后，其余的列标依次作为datExpr0的行标，但本步可以不用，因为在定义datExpr0时就已经完成了这不
datExpr0[1:6,1:6]



#########检查过度缺失值和离群样本#########

# 检查缺失值太多的基因和样本
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK




if(!gsg$allOK)
{
  #(可选)打印被删除的基因和样本名称:
  if(sum(!gsg$goodGenes)>0)
    printFlush(paste("Removinggenes:",paste(names(datExpr0)[!gsg$goodGenes], collapse =",")));
  if(sum(!gsg$goodSamples)>0)
    printFlush(paste("Removingsamples:",paste(rownames(datExpr0)[!gsg$goodSamples], collapse =",")));
  #删除不满足要求的基因和样本:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}





########聚类做离群样本检测##########


sampleTree = hclust(dist(datExpr0), method ="average");
# dist()表示转为数值，method表示距离的计算方式，其他种类的详见百度
sizeGrWindow(12,9)
# 绘制样本树:打开一个尺寸为12 * 9英寸的图形输出窗口
# 可对窗口大小进行调整
# 如要保存可运行下面语句
# pdf(file="Plots/sampleClustering.pdf",width=12,height=9);
par(cex = 0.6)	# 控制图片中文字和点的大小
par(mar =c(0,4,2,0))	# 设置图形的边界，下，左，上，右的页边距
plot(sampleTree, main ="Sample clustering to detectoutliers",sub="", xlab="", cex.lab = 1.5,
     cex.axis= 1.5, cex.main = 2)
# 参数依次表示：sampleTree聚类树，图名，副标题颜色，坐标轴标签颜色，坐标轴刻度文字颜色，标题颜色
# 其实只要包括sampleTree和图名即可




#可见有一个异常值。可以手动或使用自动方法删除它。选择一个**高度(height)**进行切割将删除异常样本，比如15(Fig1b)，并在该高度使用分支切割

# 绘制阈值切割线
abline(h = 100,col="red"); # 高度15，颜色红色
# 确定阈值线下的集群
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
# 以高度15切割，要求留下的最少为10个
table(clust)	# 查看切割后形成的集合
# clust1包含想要留下的样本.
keepSamples = (clust==1)	# 将clust序号为1的放入keepSamples
datExpr = datExpr0[keepSamples, ]		
# 将树中内容放入矩阵datExpr中，因为树中剩余矩阵不能直接作为矩阵处理
nGenes =ncol(datExpr)	# ncol，crow分别表示提取矩阵的列数和行数
nSamples =nrow(datExpr)


#########载入临床特征数据############

traitData = as.data.frame(trait2)
dim(traitData)	# 看看形状
names(traitData)	# 看看列标
# 删除不必要的列.
allTraits = traitData	# 将去掉第31和16列（两个不包含数据的列）后的traitData存入allTraits中
#allTraits = allTraits[,c(2, 11:36) ] # 在allTraits中保留第2,11到36列（只取样本名的性状相关数据）
# 这时的allTraits是只包含样本名和性状相关数据的矩阵
dim(allTraits)	# 看看形状
names(allTraits)	# 看看列标
# 形成一个包含临床特征的数据框
femaleSamples =rownames(datExpr)	# femaleSamples存放存放录入基因表达量的样本名称
traitRows =match(femaleSamples, allTraits$num)	# 将表达矩阵和性状矩阵中，样本名（Mice）重复的这些样本在allTraits中的行标返回给traitRows（一个数字向量）

traitData = as.matrix(traitData)
allTraits = as.matrix(allTraits)
datTraits = allTraits[traitRows, -1]	# 在allTraits中取上步的得到的这些行（行），并删除第一列，然后组成矩阵datTraits
rownames(datTraits) = allTraits[traitRows, 1]	# 因为上一步删了第一列，所以重新赋予第一列，即这些行的样本名字
datTraits = as.data.frame(datTraits)
collectGarbage() # 释放内存

save(datExpr, datTraits, file = "WGCNA0.3-dataInput.RData")


######1####
sampleTree2 = hclust(dist(datExpr), method ="average")
# 重新聚类样本

datTraits[] <- lapply(datTraits, as.numeric)

traitColors = numbers2colors(datTraits, signed = FALSE);
# 将临床特征值转换为连续颜色:白色表示低，红色表示高，灰色表示缺失
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels =names(datTraits),
                    main ="Sample dendrogram and trait heatmap")
# 在样本聚类图的基础上，增加临床特征值热图


########自动构建网络及识别模块#########

# 设置软阈值调参范围，powers是数组，包括1，2，...10,12,14，...,20
powers =c(c(1:10),seq(from = 12, to=20,by=2))
# 网络拓扑分析
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 绘图
sizeGrWindow(9, 5)	# 图片的宽度和高度
# 1行2列排列
par(mfrow =c(1,2));	# 一页多图，一页被分为一行，两列
cex1 = 0.9;
# 无标度拓扑拟合指数与软阈值的函数(左图)，下面的会用就行
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="SoftThreshold(power)",ylab="ScaleFreeTopologyModelFit,signedR^2",type="n",
     main =paste("Scaleindependence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# 这条线对应于h的R^2截止点
abline(h=0.90,col="red")
# Mean Connectivity与软阈值的函数(右图)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="SoftThreshold(power)",ylab="MeanConnectivity", type="n",
     main =paste("Meanconnectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],labels=powers, cex=cex1,col="red")
#找第一个大于横线的power值，发现是6，那么软阈值就是6




#########一步构建网络和识别模块##########
datExpr[] <- lapply(datExpr, as.numeric)

net = blockwiseModules(datExpr,power= 9,	 # 表达矩阵，软阈值
                       TOMType ="unsigned", minModuleSize = 2,	# 数据为无符号类型，最小模块大小为30
                       reassignThreshold = 0, mergeCutHeight = 0.25,	#mergeCutHeight合并模块的阈值，越大模块越少
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase ="femaleMouseTOM",
                       verbose = 3)



# 可视化模块
sizeGrWindow(12, 9)
# 将标签转换为颜色
mergedColors = labels2colors(net$colors)
# 绘制树状图和模块颜色图
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#保存模块赋值和模块特征基因信息，以供后续分析

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
#save(MEs, moduleLabels, moduleColors, geneTree, file="FemaleLiver-02-networkConstruction-auto.RData")



geneClust = as.data.frame(moduleLabels)
geneClust$genes = rownames(geneClust)


library(openxlsx)
#write.xlsx(geneClust, file = "geneClust.xlsx", colNames = TRUE)



##########分步法完成网络构建和模块检测#########
rm(list = ls())
library('WGCNA')


#载入第一步中的表达量和表型值
lnames = load(file = "WGCNA0.3-dataInput.RData")
lnames

# 设置网络构建参数选择范围
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
powers
# 计算无尺度分布拓扑矩阵
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#pdf(" Analysis of network topology for various soft-thresholding powers.pdf",width = 9, height=5)
#一页多图，一行2列
par(mfrow = c(1,2))
#字号
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power（左图）
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red") 
#查看最佳软阈值
sft$powerEstimate
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power（右图）
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#共表达相似性和邻接性
#通过上述计算得知，最佳软阈值为6
softPower = 6
adjacency = adjacency(datExpr, power = softPower)


#拓扑覆盖矩阵TOM
#为了减小噪声和虚假关联的影响，我们将邻接变换为拓扑重叠矩阵，并计算相应的不相似度。
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

#通过TOM进行聚类
#现在使用层次聚类来产生基因的层次聚类树(树状图)。
#注意，我们使用了hclust函数，它提供了比标准hclust函数更快的分层集群例程。


# 调用层次聚类函数
geneTree = hclust(as.dist(dissTOM), method = "average")
# 聚类结果绘图
#pdf(" Gene clustering on TOM-based dissimilarity.pdf",width = 12, height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

#在图中，每个短分枝代表一个基因，相邻的基因表示他们高度共表达。聚类方法是dynamicTreeCut包中的Dynamic TreeCut。

# 由于基因数目比较多，我们喜欢大的模块，所以我们把最小模块的大小设置得比较高，这个数值要根据自己的数据进行修改
minModuleSize = 9
# 模块识别使用动态树切割，这里敏感度deepSplit我设置了1，需要根据自己数据进行调试，默认是2
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 1, pamRespectsDendro = FALSE, minClusterSize = minModuleSize) 
# 查看模块数
table(dynamicMods)
# 转换数字标签为颜色
dynamicColors = labels2colors(dynamicMods)
# 查看颜色
table(dynamicColors)
# 绘图
#pdf("Gene dendrogram and module colors.pdf",width = 8, height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

dev.off()

#合并表达模式相似的模块
#量化各个模块的共表达相似性，根据模块间的相关性合并表达模式相似的模块。

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#  Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#  Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# 画图
#pdf("Clustering of module eigengenes.pdf",width = 7, height=6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")


#这里选择切割高度0.3，即模块之间相关性达到0.7进行合并，这个系数也要根据自己的数据进行更改，常见的是0.25（Fig.4）。

MEDissThres = 0.7
# 在树状图中加入切割线
abline(h=MEDissThres, col = "red")
# 调用自动归并函数
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# 合并模块的颜色
mergedColors = merge$colors
# 合并后新模块
mergedMEs = merge$newMEs

#为了查看合并对模块颜色的影响，我们再次绘制了基因树状图，下面依次是原始和合并后的模块颜色（Fig.5）。

#pdf("PlotsgeneDendro3.pdf",width = 12, height=9)
#pdf("Figure2_B.pdf",width = 8, height=4)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#在随后的分析中，我们将在mergedColors中使用合并的模块颜色，保存相关的变量，以便在后续分析中使用。

# 对模块颜色重命名
moduleColors = mergedColors
# 将颜色转换为数据标签
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# 保存模块颜色和标签，以供后续部分使用
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-stepByStep.RData")


########基因模块与性状关联识别重要基因#######
rm(list = ls())
library('WGCNA')
# 允许R语言以最大线程运行
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 载入第一步中的表达量和表型值
lnames = load(file = "WGCNA0.3-dataInput.RData")
lnames
# 载入第二步的网络数据
lnames = load(file = "networkConstruction-stepByStep.RData")
lnames

#将基因模块与表型数据关联
  #量化基因模块与表型数据的相关关系
  #在此分析中，我们希望识别与表型显著相关的基因模块，因为前期已经有了每个模块的概要文件(特征基因)，
  #我们只需要简单地将特征基因与表型联系起来，并寻找最显著的关联即可。




# 明确基因和样本数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# 用颜色标签重新计算MEs
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 查看每个模块的颜色及包含的基因数目
table(moduleColors)




#绘制基因模块和性状的相关性热图（Fig.1）

#pdf("Module-trait associations.pdf",width = 8, height=10)
#pdf("Figure2_C.pdf",width = 8, height=4)
# 通过p值显示相关性
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# 通过热图显示相关性
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(550), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
# 想改变热图配色可以更改 colors = greenWhiteRed(50)
dev.off()


######基因网络图######
load('WGCNA_data.RData')

genes = genes$substanceBXH

genes_type = cbind(as.data.frame(genes) , as.data.frame(moduleColors), as.data.frame(moduleLabels))

rm(list=ls()[which(ls()!='genes_type')])



load('PPIN_SYMBOL.RData')

#PPIN_SYMBOL = subset(PPIN_SYMBOL, PPIN_SYMBOL$combined_score >= 400)

PPIN = PPIN_SYMBOL[,c(15:16)]

rm(list = c('PPIN_SYMBOL', 'proInfo','genes'))

PPIN_Out = subset(PPIN, PPIN$source %in% genes_type$genes & PPIN$target %in% genes_type$genes)
PPIN_gene = rbind(as.matrix(PPIN_Out$source), as.matrix(PPIN_Out$target))
PPIN_gene = as.matrix(PPIN_gene[!duplicated(PPIN_gene)])



rm(list = c('PPIN'))

save.image('WGCNA_ppin.rdata')


######绘制基因网络图######
rm(list = ls())

load('WGCNA_ppin.rdata')

library(tidyverse)
library(igraph)


edges = PPIN_Out


vertices = as.data.frame(genes_type)


names(edges)[1:2] = c('node1', 'node2')
names(vertices)[1:2] = c('node','type')

vertices = vertices[which(vertices$node %in% PPIN_gene),]


g <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE )
#由于本示例是个无向网络，故无出度和入度之分
V(g)$degree <- degree(g)
V(g)$degree

#查看度分布
#可观察到微生物相关网络通常服从幂律分布，这个下节再讲怎样通过计算验证
degree_dist <- degree.distribution(g)[-1]
degree_num <- 1:max(V(g)$degree)

#par(mfrow = c(1, 2))
hist(V(g)$degree, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_num, degree_dist, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-intensity', main = 'Log-log degree distribution')


#加权度（Weighted degree）
V(g)$weight_degree <- strength(g)
V(g)$weight_degree


V(g)$color <- if_else(V(g)$type == "blue", "#4DBBD5B2", if_else(V(g)$type == "brown", "#E64B35B2", if_else(V(g)$type == "grey", "#bebebe", "turquoise")))

V(g)$size <- V(g)$weight_degree %>% log()*2.3


E(g)$weight <- abs(rnorm(n = ecount(g)))
E(g)$arrow.size <- .2
E(g)$edge.color <- "#FAFAFA"
graph_attr(g, "layout") <- layout_with_lgl



### 4.1 准备网络图布局数据
#?layout_in_circle # 帮助信息中，有其它布局函数。
layout1 <- layout_in_circle(g) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(g) # fr布局。
layout3 <- layout_on_grid(g) # grid布局
layout4 <- layout_in_circle(g)
#layout5 <- layout_with_mds(g)
layout6 <- layout_on_sphere(g)
layout7 = layout_with_fr(g)
layout8 = layout_nicely(g)
layout9 = layout_with_gem(g)


#pdf("Figure2_D.pdf",family = "Times",width = 10,height = 12)
#par(mar=c(5,2,1,2))


plot.igraph(g, layout=layout8,#更多参数设置信息?plot.igraph查看。
            ##节点颜色设置参数##
            vertex.color=V(g)$point.col,
            vertex.frame.color ="black",
            vertex.border=V(g)$point.col,
            ##节点大小设置参数##
            vertex.size=V(g)$n,
            ##节点标签设置参数##
            vertex.label=g$node,
            vertex.label.cex=0.8,
            #vertex.label.dist=0, # 标签距离节点中心的距离，0表示标签在节点中心。
            #vertex.label.degree = 0, # 标签对于节点的位置，0-右，pi-左，-pi/2-上，pi/2-下。
            vertex.label.col="black",
            ##链接属性参数以edge*开头##
            edge.arrow.size=0.1,
            edge.curved = TRUE
)


##设置图例，与plot.igraph()函数一起运行##
legend(
  title = "Class",
  list(x = min(layout1[,1])-0.2,
       y = min(layout1[,2])-0.17), # 图例的位置需要根据自己的数据进行调整，后面需要使用AI手动调整。
  legend = c(unique(V(g)$type)),
  fill = c("#E64B35B2","#4DBBD5B2","turquoise","#bebebe"),
  #pch=1
)


legend(
  title = "Degree",
  list(x = min(layout1[,1])+0.4,
       y = min(layout1[,2])-0.17),
  legend = c(1,seq(0,8,2)[-1]),# max(V(g)$degree)
  col = "black",
  pch=1,
  pt.lwd=1,
  pt.cex=c(1,seq(0,8,2)[-1])
)

#dev.off()

#par(mfrow = c(1, 1))


#write.csv(edges, 'genes.csv', quote = F, row.names = F, fileEncoding = 'GBK')

rm(list=ls()[which(ls()!='vertices')])
save.image('WGCNA_genes.rdata')


######富集######
rm(list = ls())

load('WGCNA_genes.rdata')


options(warn= -1)
#warnings('off')
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(enrichplot)
library(DOSE)
library(openxlsx)
library(UpSetR)
library(flexclust)
library(proxy)
library(ggplot2)
library(factoextra)


#gene = read.csv('stapAshao.csv',header = T)
gene = as.data.frame(vertices[,1])
#gene2EID = gene
#gene = as.data.frame(gene$all)
SYM2EID = function(num){
  lie = as.character(gene[,num])
  keytypes(org.Hs.eg.db)
  gene.df <- bitr(lie, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                  toType = c("ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                  OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
  #gene2EID[,num] = as.character(gene.df[,2])
  head(gene.df)
  return(gene.df)
}


######go单列#####
go = as.data.frame(c())
catch = SYM2EID(1)[,2]

ego = enrichGO(
  gene = catch,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE,
  pool = FALSE
)
#ego = DOSE::setReadable(ekegg,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
catch = as.data.frame(ego)

catch2 = catch[which(catch$ONTOLOGY == 'BP'),]
catch2 = catch2[which(catch2$Count > (fivenum(catch2$Count)[4])),]
catch2 = catch2[which(catch2$Count > (fivenum(catch2$Count)[4])),]
catch2 = catch2[which(catch2$Count > (fivenum(catch2$Count)[4])),]
catch2 = catch2[which(catch2$Count >= 16),]

catch3 = catch[which(catch$ONTOLOGY == 'MF'),]
catch3 = catch3[which(catch3$Count > (fivenum(catch3$Count)[4])),]
catch3 = catch3[which(catch3$Count >= 10),]


catch4 = catch[which(catch$ONTOLOGY == 'CC'),]
catch4 = catch4[which(catch4$Count > (fivenum(catch4$Count)[4])),]
catch4 = catch4[which(catch4$Count >= 10),]

catch_out = rbind(catch2,catch3,catch4)

library(openxlsx)
#write.xlsx(catch_out, file = "WGCNA_go.xlsx", colNames = TRUE)

#####和弦图######
rm(list = ls())
library(circlize)
library(statnet)

go = read.csv('go2.csv',header = T)

#行转列
a = go
lei = as.numeric(ncol(a))
hang = as.numeric(nrow(a))

b= c()
b = as.data.frame(b)


zs = as.numeric(lei*hang)



for (i in 1:hang) {
  for (ii in 1:(lei-1)) {
    xu = as.numeric(as.numeric(nrow(b)) + 1)
    b[xu,1] = a[i,1]
    b[xu,2] = a[i,ii+1]
  }
  print(i)
}

go_node = b[which(b$V2 != ''),]

go_node$value = 1
names(go_node)[1:2] = c('from','to')

color=NULL
#color[c("A","B","C","D","E")]=c("blue","red","yellow","green","pink")
color[go$Description]=c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")


pdf(file="Figure3_C.pdf", width=9, height=3, pointsize=8)

# 设置角度
circos.par(start.degree = 90)



chordDiagram(go_node, 
             grid.col =color,#颜色
             annotationTrack = "grid",
             transparency = 0.2,#透明度
             link.lwd = 0.00001,#线条宽度
             link.lty = 1,    # 线路类型
             link.border = 0,#边框颜色
             directional = -1,#表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
             diffHeight = mm_h(3),#外圈和中间连线的间隔
             direction.type = c("diffHeight","arrows"), #线条是否带有箭头
             link.arr.type = "big.arrow",#箭头类型
             annotationTrackHeight = c(0.04, 0.1))#网格高度

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = T, adj = c(-0.5, 0.5), cex = 0.8)
}, bg.border = NA)

goT = read.csv('go.csv',header = T)

legend("right",pch=20,legend=goT$Description,
       col=color[go$Description],bty="n",
       cex=1,pt.cex=3,border="black")

dev.off()
#####KEGG单列#########
kegg = as.data.frame(c())
catch = SYM2EID(1)[,2]

ekegg = enrichKEGG(
  gene = catch,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
ekegg = DOSE::setReadable(ekegg,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
catch = as.data.frame(ekegg)
kegg = catch
library(openxlsx)
#write.xlsx(kegg, file = "kegg.xlsx", colNames = TRUE)

rm(list=ls()[which(ls()!='kegg' &ls()!='vertices')])
#save.image('WGCNA_kegg.rdata')



library(ggplot2)
library(ggtech)
library(ggsci)



load('WGCNA_kegg.rdata')


kegg_index = read.csv('kegg_index.csv',header = F)

kegg = kegg[which(kegg$ID %in% kegg_index[,1]),]


x= kegg$p.adjust
y = factor(kegg$Description, levels = kegg$Description)

pr = ggplot(kegg, aes(x,y)) + geom_point() + geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="#ADC2DE",high ="#F78C7A") + 
  labs(x="P.adj-Value",y="Pathway name",title="KEGG Pathway enrichment")
pr







######桑葚图######

rm(list = ls())
load('Gene_MM.Rdata')

kegg = read.csv('kegg.csv',header = T)
kegg_index = read.csv('kegg_index.csv',header = F)

kegg = kegg[which(kegg$ID %in% kegg_index[,1]),]

kegg = kegg[,-1]
#行转列
a = kegg
lei = as.numeric(ncol(a))
hang = as.numeric(nrow(a))

b= c()
b = as.data.frame(b)


zs = as.numeric(lei*hang)



for (i in 1:hang) {
  for (ii in 1:(lei-1)) {
    xu = as.numeric(as.numeric(nrow(b)) + 1)
    b[xu,1] = a[i,1]
    b[xu,2] = a[i,ii+1]
  }
  print(i)
}

kegg_node = b[which(b$V2 != ''),]





library(networkD3)
library(tidyverse)


names(Gene_MM) = c('V1','V2')



phone = rbind(Gene_MM,kegg_node)
phone$V3 = 1
names(phone) = c("source","destination","n.call")


load('WGCNA_ppin.rdata')

names(vertices)[1] = c('source')
vertices = vertices[,c(1,2)]
phone = merge(phone, vertices,all.x = T)
phone$source[which(phone$type != 'NA')] = phone$type[which(phone$type != 'NA')]
phone = phone[,c(1:3)]


names(vertices)[1] = c('destination')
vertices = vertices[,c(1,2)]
phone = merge(phone, vertices,all.x = T)

phone$destination[which(phone$type != 'NA')] = phone$type[which(phone$type != 'NA')]
phone = phone[,c(1:3)]

load('MM_pinshu.RData')
MM_pinshu = MM_pinshu[which(MM_pinshu$Freq >= 343),]


shai = rbind(as.matrix(MM_pinshu$Var1),as.matrix(vertices[,1]),as.matrix(kegg$Description))
shai = as.data.frame(shai)

phone = phone[which(phone$destination %in% shai[,1] | phone$source %in% shai[,1]),]



phone = phone[!duplicated(phone),]

phone = phone[-which(phone$source == 'SNORD115@' | phone$source == 'SERPINA3'),]


phone = phone[-which(phone$destination == 'Glossitis'),]
catch=  phone[1:3,]
catch[,1] = 'Chest Pain'
catch[1,2] = 'blue'
catch[2,2] = 'turquoise'
catch[3,2] = 'grey'

phone = rbind(phone,catch)

catch = as.data.frame(table(phone$destination))
catch2 = as.data.frame(table(phone$source))

catch3 = as.data.frame(c())

catch3[1:2,1] = 'Others'
catch3[,2] = 'blue'

catch4 = as.data.frame(c())

catch4[1:2,1] = 'Others'
catch4[,2] = 'brown'
  
catch5 = as.data.frame(c())

catch5[1:9,1] = 'Others'
catch5[,2] = 'turquoise'
    
catch6 = as.data.frame(c())

catch6[1:7,1] = 'Others'
catch6[,2] = 'grey'
    
catch = rbind(catch3,catch4,catch5,catch6)
catch$n.call = 1
names(catch) = c("destination", "source","n.call" )
phone = rbind(phone,catch)





#  来源地去重
sources <- phone %>%
  distinct(source) %>%
  rename(label = source)

#目的地去重
destinations <- phone %>%
  distinct(destination) %>%
  rename(label = destination)

## 合并数据并添加一列索引
nodes <- full_join(sources, destinations, by = "label") 
nodes <- nodes %>%
  mutate(id = 1:nrow(nodes)) %>%
  select(id, everything())
head(nodes, 3)


# Rename the n.call column to weight
phone <- phone %>%
  rename(weight = n.call)

# (a) Join nodes id for source column
edges <- phone %>% 
  left_join(nodes, by = c("source" = "label")) %>% 
  rename(from = id)

# (b) Join nodes id for destination column
edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)

# (c) Select/keep only the columns from and to
edges <- select(edges, from, to, weight)
head(edges, 3)


library(tidyverse)



##注意：该包中节点和边缘列表中的ID应该是从0开始的数字值。我们直接减1 即可
nodes_d3 <- mutate(nodes, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)


library(networkD3)

#桑葚图
p = sankeyNetwork(
  Links = edges_d3, Nodes = nodes_d3, 
  Source = "from", Target = "to", 
  NodeID = "label", Value = "weight", 
  fontSize = 16, unit = "Letter(s)")

p
library(webshot)


saveNetwork(p,"sankey.html")

webshot("sankey.html" , "Figure3_A.pdf")


#write.csv(vertices, 'gene_group.csv', quote = F, row.names = F, fileEncoding = 'GBK')

