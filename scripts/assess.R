rm(list=ls())

library(bnlearn)
library(igraph)
library(graph)
type="is.mb"




load(paste("./data/trainD2C.500",type,"RData",sep="."))
load("./data/testDAG.500.RData")


BER.D2C<-NULL
BER.IAMB<-NULL
BER.GS<-NULL
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
  browser()
  cat("Dim test dataset"=dim(observedData),"\n")

  ## inference of networks with bnlearn package

  ## Ahat.GS<-(amat(gs(data.frame(observedData))))
  Ahat.IAMB<-(amat(iamb(data.frame(observedData),alpha=0.01)))
  Ahat.PC<-(amat(si.hiton.pc(data.frame(observedData),alpha=0.01)))
  Ahat.GS<-Ahat.IAMB
  igraph.GS<-graph.adjacency(Ahat.GS)
  igraph.IAMB<-graph.adjacency(Ahat.IAMB)
  igraph.PC<-graph.adjacency(Ahat.PC)


  graphTRUE<- as.adjMAT(trueDAG)
  igraph.TRUE<-graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),as.character(1:NCOL(graphTRUE))])

  ## selection of a balanced subset of edges for the assessment
  Nodes=nodes(trueDAG)
  max.edges<-min(30,length(edgeList(trueDAG)))
  subset.edges = matrix(unlist(sample(edgeList(trueDAG),size = max.edges,replace = F)),ncol=2,byrow = TRUE)
  subset.edges = rbind(subset.edges,t(replicate(n =max.edges ,sample(Nodes,size=2,replace = FALSE))))


  for(jj in 1:NROW(subset.edges)){
    i=subset.edges[jj,1]
    j=subset.edges[jj,2]
    I =as(subset.edges[jj,1],"numeric");
    J =as(subset.edges[jj,2],"numeric") ;
    pred.D2C = predict(trainD2C,I,J, observedData)

    Yhat.D2C<-c(Yhat.D2C,as.numeric(pred.D2C$response)  -1)

    Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,i,j))
    Yhat.GS<-c(Yhat.GS,is.what(igraph.GS,i,j))
    Yhat.PC<-c(Yhat.PC,is.what(igraph.PC,i,j))
    Ytrue<-c(Ytrue,is.what(igraph.TRUE,i,j)) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])

    cat(".")
  }


  ## computation of Balanced Error Rate
  BER.D2C<-BER(Ytrue,Yhat.D2C)
  BER.IAMB<-BER(Ytrue,Yhat.IAMB)
  BER.GS<-BER(Ytrue,Yhat.GS)
  BER.PC<-BER(Ytrue,Yhat.PC)
  cat("\n nDAG=",r," BER.D2C=",mean(BER.D2C), "BER.IAMB=",mean(BER.IAMB),
      "BER.GS=",mean(BER.GS),"BER.PC=",mean(BER.PC),"#0=",length(which(Ytrue==0))/length(Ytrue),"\n")

}

