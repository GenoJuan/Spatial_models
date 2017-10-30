#
# This is a script to generate random networks using BoolNet.
# Preguntas:
# ¿Qué pasa con los atractores ciclicos y complejos en la versión continua?

library(BoolNet)
library(igraph)

is.connectedRandomNetwork = function(net){
  adjMatrix = matrix(data = 0, nrow = length(net$genes), ncol = length(net$genes))
  for(i in 1:length(net$genes)){
    for(j in 1:length(net$interactions[[i]]$input)){
      adjMatrix[i,net$interactions[[i]]$input[j]] = 1
    }
  }
  adjMatrix = graph_from_adjacency_matrix(adjMatrix)
  return(is.connected(adjMatrix))
}

deFix = function(net){
  todeFix = which(net$fixed != -1)
  for(node in 1:length(todeFix)){
    net$interactions[[todeFix[node]]]$input = todeFix[node]
    net$interactions[[todeFix[node]]]$func = c(0,1)
    net$interactions[[todeFix[node]]]$expression = paste("(Gene", todeFix[node], ")", sep = "")
  }
  
  return(net)
}

total = 1e2
count = 0

while(count <= total){
  net = generateRandomNKNetwork(n = 20,
                                k = 2,
                                topology = "homogeneous",
                                #functionGeneration = 
                                noIrrelevantGenes = T,
                                readableFunctions = "short",
                                linkage = "lattice",
                                d_lattice = 1
  )
  
  plotNetworkWiring(net)
  count = count + is.connectedRandomNetwork(net)
  print(count)
  
}

if(sum(net$fixed == -1) != length(net$genes)){
  net = deFix(net = net)
}
net
attr = getAttractors(net, type = "asynchronous")
attr
