## I still keep load these two package for stand alone shiny apps

rm(list=ls())
library(D2C)
library(graph)
library(igraph)
library(readxl)
library(gRbase)
library(ROCR)
library(RBGL)



is.what<-function(iDAG,i,j,type){
  if (type=="is.mb")
    return(as.numeric(is.mb(iDAG,i,j)))
  if (type=="is.parent")
    return(as.numeric(is.parent(iDAG,i,j)))
  
  if (type=="is.child")
    return(as.numeric(is.child(iDAG,i,j)))
  if (type=="is.descendant")
    return(as.numeric(is.descendant(iDAG,i,j)))
  
  if (type=="is.ancestor")
    return(as.numeric(is.ancestor(iDAG,i,j)))
}

type="is.parent"
trainD2C<-NULL
knocked<-NULL
nNode=10
maxPar=3
wgt = 0.9
nSamples=150
additive=TRUE

g<-random_dag(1:nNode,maxpar=min(nNode,maxPar),wgt)
cnt<-2

while (sum(unlist(lapply(graph::edges(g),length)))<nNode & cnt<100){
  g<-random_dag(1:nNode,maxpar =min(nNode,maxPar),wgt)
  cnt<-cnt+1
  
}
G<-graph.adjacency(as(g,"matrix"))
paste("nEdges=",length(E(G)))

if (runif(1)<0.5){
  H = function() return(H_Rn(2)) #function() return(H_sigmoid(1))
} else {
  H = function() return(H_Rn(1))
}

DAG = new("DAG.network",
          network=as_graphnel(G),H=H,additive=additive,
          weights=c(0.5,1),sdn=runif(1,0.2,0.5))

observationsDAG <- compute(DAG,N=nSamples)#,knocked)
print(dim(observationsDAG))

# Adjust vertex size according user input
V(G)$size = 30

# Adjust arrow size according user input
E(G)$arrow.size = 0.75

# To avoid plot without boundary error
if(vcount(G) > 0){
  plot(G)
}

file="trainD2Cp.5000.20.is.parent.RData"
filetrain=paste("./data/",file,sep="")
resultsfile<-paste("./data/","ass.",file,sep="")

load(filetrain)

cat("D2C: \n # descriptors=",NCOL(allD2C@origX), 
      "\n # samples=",NROW(allD2C@origX), "\n # positives=",
    length(which(allD2C@Y==1)) , "\n")

DAG=as_graphnel(G)
Nodes=nodes(DAG)
max.edges<-length(edgeList(DAG))

if (type=="is.parent"){
  subset.edges = matrix(unlist(edgeList(DAG)),ncol=2,byrow = TRUE)
  subset.edges = unique(rbind(subset.edges,t(replicate(n =3*max.edges ,
                                                       sample(Nodes,size=2,replace = FALSE)))))
} else {
  subset.edges = unique(t(replicate(n =4*max.edges ,sample(Nodes,size=2,replace = FALSE))))
}

Yhat.D2C<-NULL
phat.D2C<-NULL
Ytrue<-NULL

for(jj in 1:NROW(subset.edges)){
  i=subset.edges[jj,1]
  j=subset.edges[jj,2]
  I =as(subset.edges[jj,1],"numeric")
  J =as(subset.edges[jj,2],"numeric") 
  if (length(intersect(c(I,J),knocked))==0){
    pred.D2C.rr =predict(trainD2C,I,J, observationsDAG,rep=4)$prob
    Yhat.D2C<-c(Yhat.D2C,round(pred.D2C.rr))
    phat.D2C<-c(phat.D2C,pred.D2C.rr)
    Ytrue<-c(Ytrue,is.what(G,i,j,type)) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
    cat(".")
  }
}
cat("\n BER=",round(BER(Ytrue,Yhat.D2C),2),"\n")
cat("AUC=",round(AUC(Ytrue,phat.D2C),2),"\n")

A=table(Ytrue,round(Yhat.D2C))
rownames(A)=c("N","P")
colnames(A)=c("N'","P'")
print(A)



