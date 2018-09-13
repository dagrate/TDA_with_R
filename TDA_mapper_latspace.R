library(TDA)
library(TDAmapper)
library(ggplot2)
library(igraph)
library(locfit)
library(ks)
library(networkD3)


# 1 <- distance function
# 2 <- distance to measure
# 3 <- kernel distance 
# 4 <- kernel density estimator (prefered)
# 5 <- k-Nearest Neighbor (last choice)
FILTER <- 4
NUMINTERVALS <- 7
PERCENT_OVERLAP <- 80
NUM_BINS_WHEN_CLUSTERING <- 20
REPULSION <- 100
BOTTLENECKDIST <- 1




# open the file
IDICT <- 2
dict <- c("VAE_normal_samplesoriginal.csv", "VAE_normal_samplesgenerated.csv")
#dict <- c("WAE_normal_samplesoriginal.csv", "WAE_normal_samplesgenerated.csv")
flnm <- dict[IDICT]
xfrauds <- read.table( paste("/home/jeremy/Documents/SnT/Data/credit_card_fraud/", flnm, sep="") , sep=',')
if(dim(xfrauds)[1]>800) {
  xfrauds <- xfrauds[1:800,]
}


lastcolumn <- dim(xfrauds)[2] - 1
colnames(xfrauds)[lastcolumn+1] <- "cc"
xfrauds[,1:lastcolumn] <- scale(xfrauds[,1:lastcolumn],center=FALSE)

xfrauds.dist = dist(xfrauds[,1:lastcolumn])


# apply kde filter
if (FILTER == 1) {
  # distance function
  xfrauds.filter <- distFct(xfrauds[,1:lastcolumn], 
                        diag(1,nrow=lastcolumn))
} else if (FILTER == 2) {
  # distance to measure
  xfrauds.filter <- dtm(xfrauds[,1:lastcolumn], 
                    diag(1,nrow=lastcolumn), 0.1)
} else if (FILTER == 3) {
  # kernel distance 
  xfrauds.filter <- kernelDist(xfrauds[,1:lastcolumn], 
                           diag(1,nrow=lastcolumn), 0.1)
} else if (FILTER == 4) {
  # kernel density estimator
  xfrauds.filter <- kde(xfrauds[,1:lastcolumn], 
                    H=diag(1,nrow=lastcolumn), 
                    eval.points = xfrauds[,1:lastcolumn])$estimate
  
} else if (FILTER == 5) {
  # k-Nearest Neighbor
  xfrauds.filter <- knnDE(xfrauds[,1:lastcolumn], 
                      diag(1,nrow=lastcolumn), 1)
}


# apply initial mapper
xfrauds.mapper <- mapper(dist_object = xfrauds.dist,
                    filter_values = xfrauds.filter,
                    num_intervals = NUMINTERVALS,
                    percent_overlap = PERCENT_OVERLAP,
                    num_bins_when_clustering = NUM_BINS_WHEN_CLUSTERING)


xfrauds.graph <- graph.adjacency(xfrauds.mapper$adjacency, mode="undirected")
#plot(xfrauds.graph )


l = length(V(xfrauds.graph))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

cc.maj.vertex <- c()
xfrauds.filter.vertex <- c()
for (i in 1:l){
  points.in.vertex <- xfrauds.mapper$points_in_vertex[[i]]
  Mode.in.vertex <- Mode(xfrauds$cc[points.in.vertex])
  cc.maj.vertex <- c(cc.maj.vertex,as.character(Mode.in.vertex))
  xfrauds.filter.vertex <- c(xfrauds.filter.vertex,mean(xfrauds.filter[points.in.vertex]))
}

vertex.size <- rep(0,l)
for (i in 1:l){
  points.in.vertex <- xfrauds.mapper$points_in_vertex[[i]]
  vertex.size[i] <- length((xfrauds.mapper$points_in_vertex[[i]]))
}

MapperNodes <- mapperVertices(xfrauds.mapper, 1:nrow(xfrauds) )
MapperNodes$cc.maj.vertex <- as.factor(cc.maj.vertex)
MapperNodes$xfrauds.filter <- xfrauds.filter.vertex
MapperNodes$Nodesize <- vertex.size

MapperLinks <- mapperEdges(xfrauds.mapper)


forceNetwork(
  Nodes = MapperNodes, 
  Links = MapperLinks, 
  Source = "Linksource", 
  Target = "Linktarget",
  Value = "Linkvalue", 
  NodeID = "Nodename",
  Group = "cc.maj.vertex", 
  #colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
  #colourScale = JS('d3.scaleOrdinal().domain(["Fraud", "Normal"]).range(["#000000", "#cccccc"])'),
  colourScale = JS('d3.scaleOrdinal().domain(["Normal"]).range(["#cccccc"])'), 
  opacity = 1, 
  linkDistance = JS("function(d){return d.value * 10}"), 
  charge = -REPULSION, # strength of the node repulsion (<0) or attraction (>0)
  legend = TRUE,
  Nodesize = "Nodesize",
  bounded = TRUE
)  



#### Rips Diagram and Barcode
if (BOTTLENECKDIST == 1) {
  rowlimit <- 800
  maxscale <- 5 # limit of the filtration
  maxdimension <- 1 # components and loops
  
  IDICT <- 1
  flnm <- dict[IDICT]
  xfrauds <- read.table( paste("/home/jeremy/Documents/SnT/Data/credit_card_fraud/", flnm, sep="") , sep=',')
  if (dim(xfrauds)[1]>rowlimit) {
    xfrauds <- xfrauds[1:rowlimit,]  
  }
  #Diag1 <- ripsDiag(X = xfrauds[,1:lastcolumn], maxdimension, maxscale, library = "GUDHI", printProgress = TRUE)
  Diag1 <- alphaComplexDiag(X = xfrauds[,1:lastcolumn], printProgress = TRUE)
  #Diag1 <- alphaShapeDiag(X = xfrauds[,1:lastcolumn], printProgress = TRUE)
  #plot(Diag1[["diagram"]], barcode = TRUE)
  
  IDICT <- 2
  flnm <- dict[IDICT]
  xfrauds <- read.table( paste("/home/jeremy/Documents/SnT/Data/credit_card_fraud/", flnm, sep="") , sep=',')
  if (dim(xfrauds)[1]>rowlimit) {
    xfrauds <- xfrauds[1:rowlimit,]
  }
  #Diag2 <- ripsDiag(X = xfrauds[,1:lastcolumn], maxdimension, maxscale, library = "GUDHI", printProgress = TRUE)
  Diag2 <- alphaComplexDiag(X = xfrauds[,1:lastcolumn], printProgress = TRUE)
  #Diag2 <- alphaShapeDiag(X = xfrauds[,1:lastcolumn], printProgress = TRUE)
  #plot(Diag2[["diagram"]], barcode = TRUE)
  
  print(bottleneck(Diag1[["diagram"]], Diag2[["diagram"]]))
}