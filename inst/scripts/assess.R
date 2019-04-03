### Script comparing inference accuracy of D2C and some "bnlearn" causal inference algorithms  (IAMB and PC)

rm(list=ls())

library(bnlearn)
library(igraph)
library(graph)
library(gRbase)

type="is.parent"
## kind of causal relationship  assessed

load(paste("./data/trainD2C.1000",type,"RData",sep="."))
load("./data/testDAG.500.RData")


BER.D2C<-NULL
BER.IAMB<-NULL

BER.PC<-NULL
Yhat.D2C<-NULL
Yhat.IAMB<-NULL
Yhat.GS<-NULL
Yhat.PC<-NULL
Ytrue<-NULL
for ( r in 1:testDAG@NDAG){
  set.seed(r)
  observedData<-testDAG@list.observationsDAGs[[r]]
  trueDAG<-testDAG@list.DAGs[[r]]
  
  cat("Dim test dataset"=dim(observedData),"\n")
  
  ## inference of networks with bnlearn package
  Ahat.IAMB<-(amat(iamb(data.frame(observedData),alpha=0.01)))
  Ahat.PC<-(amat(si.hiton.pc(data.frame(observedData),alpha=0.01)))
 
  colnames(Ahat.IAMB)<-colnames(observedData)
  
  colnames(Ahat.PC)<-colnames(observedData)
 
  igraph.IAMB<-graph.adjacency(Ahat.IAMB)
  igraph.PC<-graph.adjacency(Ahat.PC)
  
  
  graphTRUE<- as.adjMAT(trueDAG)
  igraph.TRUE<-graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),as.character(1:NCOL(graphTRUE))])
  
  ## selection of a balanced subset of edges for the assessment
  Nodes=nodes(trueDAG)
  max.edges<-min(30,length(edgeList(trueDAG)))
  if (type=="is.parent"){
    subset.edges = matrix(unlist(sample(edgeList(trueDAG),size = max.edges,replace = F)),ncol=2,byrow = TRUE)
    subset.edges = rbind(subset.edges,t(replicate(n =max.edges ,sample(Nodes,size=2,replace = FALSE))))
  } else {
    
    subset.edges = t(replicate(n =3*max.edges ,sample(Nodes,size=2,replace = FALSE)))
  }
  
  for(jj in 1:NROW(subset.edges)){
    i=subset.edges[jj,1]
    j=subset.edges[jj,2]
    I =as(subset.edges[jj,1],"numeric");
    J =as(subset.edges[jj,2],"numeric") ;
    pred.D2C = predict(trainD2C,I,J, observedData)
    
    Yhat.D2C<-c(Yhat.D2C,as.numeric(pred.D2C$response)  -1)
    
    Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,i,j,type=type))
    
    Yhat.PC<-c(Yhat.PC,is.what(igraph.PC,i,j,type=type))
    Ytrue<-c(Ytrue,is.what(igraph.TRUE,i,j,type=type)) 
    
    cat(".")
  }
  
  
  ## computation of Balanced Error Rate
  BER.D2C<-BER(Ytrue,Yhat.D2C)
  BER.IAMB<-BER(Ytrue,Yhat.IAMB)
  
  BER.PC<-BER(Ytrue,Yhat.PC)
  cat("\n nDAG=",r," BER.D2C=",mean(BER.D2C), "BER.IAMB=",mean(BER.IAMB),
      "BER.PC=",mean(BER.PC),"#0=",length(which(Ytrue==0))/length(Ytrue),"\n")
  
}

