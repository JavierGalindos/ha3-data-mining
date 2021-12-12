# Exercise 2
# Javier Galindos

rm(list = ls()) # Remove all the objects we created so far.

library("rstudioapi")  

# Set working directory
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()  


library(igraph) # Load the igraph package
library(gdata) # Distances for paths

# Loading the graphs

# DATASET 1: edgelist 
nodes <- read.csv("Datasets/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Datasets/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL


net.dir <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
net.undir <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# Remove loops and multiple
net.dir <- simplify(net.dir, remove.multiple = T, remove.loops = T)
net.undir <- simplify(net.undir, remove.multiple = T, remove.loops = T)
plot(net.dir,edge.arrow.size=.4,edge.curved=0)
plot(net.undir,edge.arrow.size=.4,edge.curved=0)


# Prepare Graph's Adjacency Matrix

graph_AM <- as_adjacency_matrix(net.undir, attr="weight")
graphMatrix <- as.matrix(graph_AM)

graph_AM.dir <- as_adjacency_matrix(net.dir, attr="weight")
graphMatrix.dir <- as.matrix(graph_AM.dir)


# Local Clustering Coefficient (Undirected graph)

JG_LCC <- function(graphMatrix){
  LCC <- matrix(,nrow = nrow(graphMatrix),1)
  nRow <- nrow(graphMatrix)
  rownames(LCC) <- rownames(graphMatrix)
  
  for(i in 1:nRow) {
    neigh <- which(graphMatrix[,i] > 0)
    
    if (length(neigh) > 0){
      nr <- length(neigh)
      count = 0
      for(j in 1 : (length(neigh) - 1)) {
        for(k in (j + 1) : length(neigh)){
          if(graphMatrix[neigh[j],neigh[k]]>0)
            count <- count+1
        }
      }
      LCC[i] <- ((2*count) / (nr*(nr-1)))
    }
    else {
      LCC[i] <- 0
    }
    
  }
  return(LCC)
}


lcc<- JG_LCC(graphMatrix)
lcc
transitivity(net.undir, type="local")

# Degree centrality (Undirected graph)

JG_DC <- function(graphMatrix) {
  DC <- matrix(0, nrow = nrow(graphMatrix), 1)
  nRow <- nrow(graphMatrix)
  rownames(DC) <- rownames(graphMatrix)
  
  for(i in 1:nRow) {
    count <- 0
    for(j in 1:nRow) {
      if(graphMatrix[i,j] > 0 & i != j) {
        count <- count+1
      }
      DC[i] <- count / (nRow - 1)
    }
  }
  return(DC)
}

dc <- JG_DC(graphMatrix)
dc
degree(net.undir, mode="total")/(length(V(net.undir))-1)


# Degree Prestige (Directed graph)

JG_DP <- function(graphMatrix) {
  DP <- matrix(0, nrow = nrow(graphMatrix), 1)
  nRow <- nrow(graphMatrix)
  rownames(DP) <- rownames(graphMatrix)
  
  for (i in 1:nRow) {
    count <- 0
    for(j in 1:nRow) {
      # In
      if(graphMatrix[j,i] > 0 & i!=j ){
        count<-count+1
      }
      DP[i] <- count / (nRow-1)
    }
  }
  return(DP)
}

dp <- JG_DP(graphMatrix.dir)
dp
degree(net.dir, mode="in")/(length(V(net.undir))-1)

# Gregariousness of a node (Directed graph)

JG_Greg <- function(graphMatrix) {
  Greg <- matrix(0, nrow = nrow(graphMatrix), 1)
  nRow <- nrow(graphMatrix)
  rownames(Greg) <- rownames(graphMatrix)
  
  for (i in 1:nRow) {
    count <- 0
    for(j in 1:nRow) {
      # Out
      if(graphMatrix[i,j] > 0 & i!=j ){
        count<-count+1
      }
      Greg[i] <- count / (nRow-1)
    }
  }
  return(Greg)
}

greg <- JG_Greg(graphMatrix.dir)
greg
degree(net.dir, mode="out")/(length(V(net.undir))-1)

# Closeness centrality (Undirected graphs) and proximity prestige (Directed graphs)

JG_CC <- function(graphMatrix){
  # Convert to a Graph object to compute distances
  Graph <- graph.adjacency(graphMatrix, weighted=NULL, mode = "undirected")
  # Get the distances of the paths
  Dist <- distances(Graph, v=V(Graph), to=V(Graph), mode = "all", weights = NULL)
  CC <- matrix(,nrow = nrow(graphMatrix), 1)
  nRow <- nrow(graphMatrix)
  rownames(CC) <- rownames(graphMatrix)
  
  for (i in 1:nRow) {
    CC[i] <- 1 / (sum(Dist[i,]) / (nRow-1))
  }
  return(CC)
}

cc <- JG_CC(graphMatrix)
cc
closeness(net.undir, mode="all", weights=NA) *(length(V(net.dir))-1)


JG_Prox <- function(graphMatrix){
  # Convert to a Graph object to compute distances
  Graph <- graph.adjacency(graphMatrix, weighted=NULL, mode = "directed")
  # Get the distances of the paths
  Dist <- distances(Graph, v=V(Graph), to=V(Graph))
  Prox <- matrix(,nrow = nrow(graphMatrix), 1)
  nRow <- nrow(graphMatrix)
  rownames(Prox) <- rownames(graphMatrix)
  
  for (i in 1:nRow) {
    # influence is defined as the follower of ith element.
    influence <- which(Dist[i,] != Inf)
    influFrac <- length(influence) / (nrow(Dist) - 1)
    avg_dist <- sum(Dist[i,influence]) / length(influence)
    Prox[i] <- influFrac / avg_dist
  }
  return(Prox)
}

prox <- JG_Prox(graphMatrix.dir)
prox



# Betweenness Centrality (Directed graphs) 

JG_BC <- function(net){
  nRow <- length(V(net))
  BC <- matrix(,nrow = nRow, 1)
  
  for (i in 1:nRow) {
    # Fraction of pairs
    f_jk <- 0
    for(j in 1:nRow) {
      for(k in 1:nRow) {
        if (j != i && k != i){
          # Shortest paths
          q_jk <- suppressWarnings(all_shortest_paths(net,j,k)$res)
          if (length(q_jk) > 0){ # Avoid division by 0
            # Paths where node i is crossed by
            vertex = i
            q_jk_i <- sum(sapply(q_jk,function(x){vertex %in% x}))
            # Fraction
            f_jk <- f_jk +  q_jk_i / length(q_jk)
          } else{
            f_jk <- f_jk + 0
          }
        }
      }
    }
    BC[i]<-f_jk / ((nRow)*(nRow-1))
  }
  return(BC)
}

bc <- JG_BC(net.dir)
bc

betweenness(net.dir, directed=T, weights=NULL)/ ((length(V(net.dir))-1)*length(V(net.dir)))



# Common neighbor based measure

JG_CN <- function(net){
  nRow <- length(V(net))
  CN <- matrix(,nRow, nRow)

  for(i in 1:nRow) {
    neigh_b <- neighbors(net,i)
    for(j in 1:nRow) {
      neigh_a <- neighbors(net,j)
      CN[i,j] <- length(intersection(neigh_a,neigh_b))
    }
  }
  return (CN)
}

cn <- JG_CN(net.undir)
cn

# Jaccard Measure
JG_Jaccard <- function(net){
  nRow <- length(V(net))
  Jaccard <- matrix(,nRow, nRow)
  
  for(i in 1:nRow) {
    neigh_b <- neighbors(net,i)
    for(j in 1:nRow) {
      neigh_a <- neighbors(net,j)
      Jaccard[i,j] <- length(intersection(neigh_a,neigh_b)) / length(union(neigh_a,neigh_b))
    }
  }
  return (Jaccard)
}

jac <- JG_Jaccard(net.undir)
jac
similarity(net.undir, method = "jaccard")
