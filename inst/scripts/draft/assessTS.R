rm(list=ls())
library(D2C)
require(bnlearn)
library(pcalg)
library(kpcalg)
library(lmtest)
type="is.parent"

is.what<-function(iDAG,i,j){
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

noNodes<-c(5,6)
## range of number of nodes

N<-c(200,300)
## range of number of samples

NDAG=5
## number of DAGs to be created and simulated
NDAG.test=10
nseries=c(20,100)
sdev<-c(0.1,0.3)

goParallel=FALSE
savefile<-FALSE
namefile<-"../D2Cdata2/traintestSTAR.200.100.RData"

load(namefile)
namefile2<-"../D2Cdata2/testSTAR.RData"
if (TRUE){
  testDAG<-new("simulatedTS",NDAG=NDAG.test, N=c(300,400), noNodes=noNodes,
               seed=101,sdn=sdev,goParallel=goParallel,nseries=nseries,
               typeser=setdiff(8:18,14))
  
  save(file=namefile2, list=c("testDAG"))
} else {
  
  load(namefile2)
}
BER.D2C<-NULL
BER.D2C.1<-NULL
BER.IAMB<-NULL
BER.GS<-NULL
BER.PC<-NULL
BER.KPC<-NULL
Yhat.D2C<-NULL
Yhat.D2C.1<-NULL
Yhat.IAMB<-NULL
Yhat.GS<-NULL
Yhat.PC<-NULL
Yhat.KPC<-NULL
Yhat.GRA<-NULL
Ytrue<-NULL
for ( r in 1:testDAG@NDAG){
  set.seed(r)
  observedData<-testDAG@list.observationsDAGs[[r]]
  N<-NROW(observedData)
  observedData<-observedData[1:min(N,300),]
  N<-NROW(observedData)
  trueDAG<-testDAG@list.DAGs[[r]]
  
  cat("Dim test dataset"=dim(observedData),"\n")
  
  ## inference of networks with bnlearn package
  
  ## Ahat.GS<-(amat(gs(data.frame(observedData))))
  Ahat.IAMB<-(amat(iamb(data.frame(observedData),optimized=TRUE,
                        alpha=0.01,max.sx=3)))
  print("Done IAMB")
  Ahat.PC<-(amat(si.hiton.pc(data.frame(observedData),alpha=0.01)))
  #Ahat.IAMB<-Ahat.PC
  print("Done Si.Hiton")
  suffStat <- list(C = cor(observedData),n=NROW(observedData))
  normal.pag <- pc(suffStat, indepTest=gaussCItest, alpha = 0.01, ,m.max=5,
                   verbose=FALSE,p=NCOL(observedData),
                   numCores=3)
  Ahat.GS<-as(normal.pag, "matrix") 
  print("Done pc")
  # kpag <- kpc(suffStat = list(data=observedData, ic.method="hsic.gamma"),
  #              indepTest=kernelCItest, alpha = 0.01, 
  #                   verbose=FALSE,p=NCOL(observedData),m.max=3)
  Ahat.KPC<-Ahat.GS ##as(kpag, "matrix")
  print("Done kpc")
  igraph.GS<-igraph::graph.adjacency(Ahat.GS)
  igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
  igraph.PC<-igraph::graph.adjacency(Ahat.PC)
  igraph.KPC<-igraph::graph.adjacency(Ahat.KPC)
  
  graphTRUE<- gRbase::as.adjMAT(trueDAG)
  igraph.TRUE<-igraph::graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),as.character(1:NCOL(graphTRUE))])
  
  ## selection of a balanced subset of edges for the assessment
  Nodes=graph::nodes(trueDAG)
  max.edges<-min(30,length(gRbase::edgeList(trueDAG)))
  subset.edges = matrix(unlist(sample(gRbase::edgeList(trueDAG),size = max.edges,replace = F)),
                        ncol=2,byrow = TRUE)
  subset.edges = rbind(subset.edges,t(replicate(n =max.edges ,
                                                sample(Nodes,size=2,replace = FALSE))))
  
  
  for(jj in 1:NROW(subset.edges)){
    i =as(subset.edges[jj,1],"numeric");
    j =as(subset.edges[jj,2],"numeric") ;
    pred.D2C = predict(trainD2C.1,i,j, observedData)
    pred.D2C.1 = pred.D2C #predict(trainD2C.1,i,j, observedData)
    pred.GRANGER=pred.D2C.1 #1-grangertest(observedData[,i],observedData[,j], order = 1)$"Pr(>F)"[2]
    Yhat.GRA<-c(Yhat.GRA,as.numeric(pred.GRANGER>0.5) )
    Yhat.D2C<-c(Yhat.D2C,as.numeric(pred.D2C$response)  -1)
    Yhat.D2C.1<-c(Yhat.D2C.1,as.numeric(pred.D2C.1$response)  -1)
    Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,i,j))
    Yhat.GS<-c(Yhat.GS,is.what(igraph.GS,i,j))
    Yhat.PC<-c(Yhat.PC,is.what(igraph.PC,i,j))
    Yhat.KPC<-c(Yhat.KPC,is.what(igraph.KPC,i,j))
    Ytrue<-c(Ytrue,is.what(igraph.TRUE,i,j)) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
    
    cat(".")
  }
  
  
  ## computation of Balanced Error Rate
  BER.D2C<-BER(Ytrue,Yhat.D2C)
  BER.D2C.1<-BER(Ytrue,Yhat.D2C.1)
  BER.IAMB<-BER(Ytrue,Yhat.IAMB)
  BER.GS<-BER(Ytrue,Yhat.GS)
  BER.PC<-BER(Ytrue,Yhat.PC)
  BER.KPC<-BER(Ytrue,Yhat.KPC)
  BER.GRA<-BER(Ytrue,Yhat.GRA)
  cat("\n r=",r," BER.D2C=",mean(BER.D2C), " BER.D2C.1=",mean(BER.D2C.1),
      "BER.IAMB=",mean(BER.IAMB),"BER.GRA=",mean(BER.GRA),
      "BER.GS=",mean(BER.GS),"BER.PC=",mean(BER.PC),"BER.KPC=",mean(BER.KPC),
      "#0=",length(which(Ytrue==0))/length(Ytrue),"\n")
  
}

