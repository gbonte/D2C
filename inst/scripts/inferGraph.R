
rm(list=ls())
library(D2C)
library(graph)
library(igraph)
library(gRbase)
library(ROCR)
library(RBGL)
library(bnlearn)


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

set.seed(0)

type="is.parent"
trainD2C<-NULL
knocked<-NULL
nNode=20
maxPar=4
wgt = 0.9
nSamples=250
additive=FALSE

g<-random_dag(1:nNode,maxpar=min(nNode,maxPar),wgt)
cnt<-0

while (sum(unlist(lapply(graph::edges(g),length)))<nNode & cnt<100){
  g<-random_dag(1:nNode,maxpar =min(nNode,maxPar),wgt)
  cnt<-cnt+1
  
}
G<-graph.adjacency(as(g,"matrix"))


if (runif(1)<0.5){
  H = function() return(H_Rn(2)) #function() return(H_sigmoid(1))
} else {
  H = function() return(H_Rn(1))
}

DAG = new("DAG.network",
          network=as_graphnel(G),H=H,additive=additive,
          weights=c(0.5,1),sdn=runif(1,0.2,0.5))

observationsDAG <- compute(DAG,N=nSamples)
cat("Observed data=",dim(observationsDAG), 
    "nEdges=",length(E(G)),"\n")

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

cat("D2C learner: \n # descriptors=",NCOL(allD2C@origX), 
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


Ahat.IAMB<-(amat(iamb(data.frame(observationsDAG),alpha=0.01,max.sx=3)))
igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)


Yhat.D2C<-NULL
Yhat.IAMB<-NULL
phat.D2C<-NULL
Ytrue<-NULL
cat("D2C inferring", NROW(subset.edges), 
    "direct dependencies: please wait \n")

for(jj in 1:NROW(subset.edges)){
  i=subset.edges[jj,1]
  j=subset.edges[jj,2]
  I =as(subset.edges[jj,1],"numeric")
  J =as(subset.edges[jj,2],"numeric") 
  if (length(intersect(c(I,J),knocked))==0){
    pred.D2C.rr =predict(trainD2C,I,J, observationsDAG,rep=4)$prob
    Yhat.D2C<-c(Yhat.D2C,round(pred.D2C.rr))
    phat.D2C<-c(phat.D2C,pred.D2C.rr)
    Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,I,J,type))
    Ytrue<-c(Ytrue,is.what(G,i,j,type)) 
    cat(".")
  }
}



cat("\n Assessment accuracy: \n BER D2C=",
    round(BER(Ytrue,Yhat.D2C),2),"\n")
cat("AUC D2C=",round(AUC(Ytrue,phat.D2C),2),"\n")

cat("\n  BER IAMB=",
    round(BER(Ytrue,Yhat.IAMB),2),"\n")

A=table(Ytrue,round(Yhat.D2C))
rownames(A)=c("N","P")
colnames(A)=c("N'","P'")
print(A)



