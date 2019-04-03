### Script comparing D2C and ParallelPC

rm(list=ls())

library(bnlearn)
library(igraph)
library(graph)
library(gRbase)
require(ROCR)
require(pcalg)
require(ParallelPC)
type="is.parent"





load(url("https://dl.dropboxusercontent.com/u/15579987/D2Cdata/train.create.D2C.10000.is.parent.RData"))
print(paste("# descriptors=",NCOL(trainD2C@X), "\n # samples=",NROW(trainD2C@X), "\n # positives=",length(which(trainD2C@Y==1)) ))
nKnockeDown=0

nSamples=150
wgt = 0.9
cnodes=0
cedges=0
Nnodes=NULL
Nedges=NULL
aAUC=NULL
for (nNode in seq(20,100,by=5)){
  
  
  for (maxPar in seq(4,10,by=2)){
    knocked<-NULL
    if (nKnockeDown>0)
      knocked<<-sample(1:nNode,min(nNode,nKnockeDown))
    
    g<-random_dag(1:nNode,maxpar=min(nNode,maxPar),wgt)
    
    cnt<-2
    
    while (sum(unlist(lapply(graph::edges(g),length)))<nNode & cnt<100){
      g<-random_dag(1:nNode,maxpar =min(nNode,maxPar),wgt)
      cnt<-cnt+1
      
    }
    G<-graph.adjacency(as(g,"matrix"))
    if (!is.null(G))
      print(paste("nNodes=",nNode, ": nEdges=",length(E(G))))
    Nnodes=c(Nnodes,nNode)
    Nedges=c(Nedges,length(E(G)))
    if (runif(1)<0.5){
      H = function() return(H_sigmoid(1))
    } else {
      H = function() return(H_Rn(1))
    }
    
    additive=sample(c(TRUE,FALSE),1)
    
    DAG = new("DAG.network",
              network=as_graphnel(G),H=H,additive=additive,sdn=runif(1,0.2,0.5))
    
    observationsDAG <<- compute(DAG,N=nSamples,knocked)
    
    p<-ncol(observationsDAG)
    suffStat<-list(C=cor(observationsDAG),n=nrow(observationsDAG))
    P<-igraph.from.graphNEL(pc_parallel(suffStat, indepTest=gaussCItest, p=p, skel.method="parallel", alpha=0.01, num.cores=2)@graph)
    
    
    DAG=as_graphnel(G)
    Nodes=nodes(DAG)
    max.edges<-length(edgeList(DAG))
    
    if (type=="is.parent"){
      subset.edges = matrix(unlist(edgeList(DAG)),ncol=2,byrow = TRUE)
      subset.edges = unique(rbind(subset.edges,t(replicate(n =max.edges ,sample(Nodes,size=2,replace = FALSE)))))
    } else {
      subset.edges = unique(t(replicate(n =4*max.edges ,sample(Nodes,size=2,replace = FALSE))))
    }
    
    Yhat.D2C<-NULL
    Yhat.PC<-NULL
    phat.D2C<-NULL
    Ytrue<-NULL
    
    for(jj in 1:NROW(subset.edges)){
      i=subset.edges[jj,1]
      j=subset.edges[jj,2]
      I =as(subset.edges[jj,1],"numeric")
      J =as(subset.edges[jj,2],"numeric") 
      if (length(intersect(c(I,J),knocked))==0){
        pred.D2C.rr<-NULL
        for (rr in 1:2){
          Irr<-sample(NROW(observationsDAG),NROW(observationsDAG)-2)
          pred.D2C.rr =c(pred.D2C.rr, predict(trainD2C,I,J, observationsDAG[Irr,])$prob[1,"1"])
        }
        Yhat.D2C<-c(Yhat.D2C,round(mean(pred.D2C.rr)))
        
        phat.D2C<-c(phat.D2C,mean(pred.D2C.rr))
        Ytrue<-c(Ytrue,is.what(G,i,j,type)) 
        Yhat.PC<-c(Yhat.PC,is.what(P,i,j,type)) 
        cat(".")
      }
      
      
    }
    
    
    cat(paste("\n"," nNodes=",nNode, ": nEdges=",length(E(G)),": BER=",round(BER(Ytrue,Yhat.D2C),2),
              ": AUC=",round(AUC(Ytrue,phat.D2C),2),"BER.parallePC=",round(BER(Ytrue,Yhat.PC),2)))
    
    aAUC=c(aAUC,round(AUC(Ytrue,phat.D2C),2))
  }
}