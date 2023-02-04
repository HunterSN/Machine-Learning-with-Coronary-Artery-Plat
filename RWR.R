#rm(list = ls())

#load('PPIN_network.RData')

#rm(list=ls()[which(ls()!='MM_CAQ_Net' & ls()!='MM_target_CAQ_Net'& ls()!='MM_target_Net'& ls()!='PPIN_MM_target_Net'& ls()!='PPIN_Out'& ls()!='PPIN_MM_target_CAQ_Net')])

#MM_CAQ_Net = as.data.frame(MM_CAQ_Net)
#MM_target_CAQ_Net = as.data.frame(MM_target_CAQ_Net)
#MM_target_Net = as.data.frame(MM_target_Net)
#PPIN_Out = as.data.frame(PPIN_Out)
#PPIN_MM_target_CAQ_Net = as.data.frame(PPIN_MM_target_CAQ_Net)


#save.image(file = "net_lianglie.RData")

rm(list = ls())

load('net_lianglie.RData')

library(RandomWalkRestartMH)
library(igraph)



#######构造igraph关系网络(简单版)#######


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





####单网络####

rm(list = ls())

load('net_lianglie.RData')

library(RandomWalkRestartMH)
library(igraph)




rm(list=ls()[which(ls()!='PPIN_MM_target_CAQ_Net')])


data = PPIN_MM_target_CAQ_Net

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


PPIN_MM_target_CAQ_Net_index = as.data.frame(cbind(as.matrix(labels), as.matrix(ids)))
names(PPIN_MM_target_CAQ_Net_index) = c('labels', 'NodeNames')
rownames(PPIN_MM_target_CAQ_Net_index) = NULL





PPI_MultiplexObject <- create.multiplex(list(PPI=g))


AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)


#save.image(file = "RWRSource_PPIN_MM_target_CAQ_Net.RData")
##########RWRSource_PPIN_MM_target_CAQ_Net###############

rm(list = ls())

load('RWRSource_PPIN_MM_target_CAQ_Net.RData')


library(RandomWalkRestartMH)
library(igraph)


seed = subset(PPIN_MM_target_CAQ_Net_index, PPIN_MM_target_CAQ_Net_index$labels == 'Coronary artery plaque')

SeedGene <- seed[1,2]
## We launch the algorithm with the default parameters (See details on manual)
RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,
                                                 PPI_MultiplexObject,SeedGene)
# We display the results
RWR_PPI_Results


## In this case we selected to induce a network with the Top 15 genes.
TopResults_PPI <-
  create.multiplexNetwork.topResults(RWR_PPI_Results,PPI_MultiplexObject,
                                     k=15)


## We print that cluster with its interactions.
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI, vertex.label.color="black",vertex.frame.color="#ffffff",
     vertex.size= 20, edge.curved=.2,
     vertex.color = ifelse(igraph::V(TopResults_PPI)$name == "PIK3R1","yellow",
                           "#00CCFF"), edge.color="blue",edge.width=0.8)



RWR_PPIN_MM_target_CAQ_seedCAQ = RWR_PPI_Results$RWRM_Results
RWR_PPIN_MM_target_CAQ_seedCAQ = merge(RWR_PPIN_MM_target_CAQ_seedCAQ, PPIN_MM_target_CAQ_Net_index, all.x = T)

rm(list=ls()[which(ls()!='RWR_PPIN_MM_target_CAQ_seedCAQ')])
#RWR_PPIN_MM_target_CAQ_seed_Myocardial_bridge = RWR_PPIN_MM_target_CAQ_seedCAQ

#rm(RWR_PPIN_MM_target_CAQ_seedCAQ)
save.image(file = "RWR_PPIN_MM_target_CAQ_seed_Myocardial_bridge.RData")

####异构网络####

rm(list = ls())

load('net_lianglie.RData')

rm(list=ls()[which(ls()!='MM_CAQ_Net' & ls()!='MM_target_Net'& ls()!='PPIN_Out')])

#MM_Net = igraph
#save(MM_Net, file = 'MM_Net.RData')
load('MM_Net.RData')



library(RandomWalkRestartMH)
library(igraph)


Disease_Network = MM_Net
GeneDiseaseRelations = MM_target_Net

names(GeneDiseaseRelations) = c("hgnc_symbol", "mim_morbid") 







#PPIN_list
data = PPIN_Out

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


PPIN_index = as.data.frame(cbind(as.matrix(labels), as.matrix(ids)))
names(PPIN_index) = c('hgnc_symbol', 'NodeNames')
rownames(PPIN_index) = NULL



PPIN_seed = subset(PPIN_index, PPIN_index$hgnc_symbol == 'ACE2')


PPI_MultiplexObject <- create.multiplex(list(PPI=g))

Disease_MultiplexObject <- create.multiplex(list(Disease=Disease_Network))



GeneDiseaseRelations = merge(GeneDiseaseRelations, PPIN_index, all.x = T)
GeneDiseaseRelations$hgnc_symbol = GeneDiseaseRelations$NodeNames
GeneDiseaseRelations = GeneDiseaseRelations[,-3]

## We keep gene-diseases associations where genes are present in the PPI
## network
GeneDiseaseRelations_PPI <-
  GeneDiseaseRelations[which(GeneDiseaseRelations$hgnc_symbol %in%
                               PPI_MultiplexObject$Pool_of_Nodes),]

## We create the MultiplexHet object.
PPI_Disease_Net <- create.multiplexHet(PPI_MultiplexObject,
                                       Disease_MultiplexObject, GeneDiseaseRelations_PPI)

PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)

SeedDisease <- c("Chest Pain")
SeedGene <- PPIN_seed[1,2]

RWRH_PPI_Disease_Results <-
  Random.Walk.Restart.MultiplexHet(PPIHetTranMatrix,
                                   PPI_Disease_Net,SeedGene,SeedDisease)


## In this case we select to induce a network with the Top 10 genes
## and the Top 10 diseases.
TopResults_PPI_Disease <-
  create.multiplexHetNetwork.topResults(RWRH_PPI_Disease_Results,
                                        PPI_Disease_Net, GeneDiseaseRelations_PPI, k=10)
## We print that cluster with its interactions.
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI_Disease, vertex.label.color="black",
     vertex.frame.color="#ffffff",
     vertex.size= 20, edge.curved=.2,
     vertex.color = ifelse(V(TopResults_PPI_Disease)$name == "PIK3R1"
                           | V(TopResults_PPI_Disease)$name == "269880","yellow",
                           ifelse(V(TopResults_PPI_Disease)$name %in% 
                                    PPI_Disease_Net$Multiplex1$Pool_of_Nodes,"#00CCFF","Grey75")),
     edge.color=ifelse(E(TopResults_PPI_Disease)$type == "PPI","blue",
                       ifelse(E(TopResults_PPI_Disease)$type == "Disease","black","grey50")),
     edge.width=0.8,
     edge.lty=ifelse(E(TopResults_PPI_Disease)$type == "bipartiteRelations",
                     2,1),
     vertex.shape= ifelse(V(TopResults_PPI_Disease)$name %in%
                            PPI_Disease_Net$Multiplex1$Pool_of_Nodes,"circle","rectangle"))




