rm(list = ls())


library(igraph)

####MM_CAQ####
load('MM&CAQ_matrix.RData')



mat = as.matrix(symMM_CAQ_matrix)
g <- graph.adjacency(mat, weighted=TRUE)


MM_CAQ_path <- shortest.paths(g, algorithm = "dijkstra")

rm(list=ls()[which(ls()!='MM_CAQ_path')])
#save.image(file = "MM_CAQ_dijkstra.RData")

#####MM_target#######
rm(list = ls())

load('RWRSource_PPIN_MM_target_CAQ_Net.RData')


library(RandomWalkRestartMH)
library(igraph)




mat = AdjMatrix_PPI
g <- graph.adjacency(mat, weighted=TRUE)


MM_CAQ_target_path <- shortest.paths(g, algorithm = "dijkstra")

rm(list=ls()[which(ls()!='MM_CAQ_path')])
#save.image(file = "MM_CAQ_dijkstra.RData")


