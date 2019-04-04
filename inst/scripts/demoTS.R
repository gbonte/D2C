rm(list=ls())
library(D2C)
require(bnlearn)
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

noNodes<-c(10,20)
## range of number of nodes

N<-c(500,1000)
## range of number of samples

NDAG=50
## number of DAGs to be created and simulated
NDAG.test=100

sdev<-c(0.1,0.3)

goParallel=FALSE
savefile<-TRUE
namefile<-paste("./data/trainD2C",NDAG,type,"RData",sep=".")
if (TRUE){
  
  
  trainDAG<-new("simulatedTS",NDAG=NDAG, N=N, noNodes=noNodes,
                seed=10,sdn=sdev,goParallel=goParallel)
  
  
  descr.example<-new("D2C.descriptor",bivariate=TRUE,ns=3,acc=TRUE,lin=FALSE)
  
  trainD2C<-new("D2C",sDAG=trainDAG,
                descr=descr.example,ratioEdges=1,
                max.features=50, type=type,goParallel=goParallel,verbose=TRUE)
  
  
  trainD2C@mod
  predict(trainD2C@mod,trainD2C@X)
  
  testDAG<-new("simulatedTS",NDAG=NDAG.test, N=N, noNodes=noNodes,
               seed=101,sdn=sdev,goParallel=goParallel)
  
  
  if (FALSE){
    testD2C<-new("D2C",sDAG=testDAG,
                 descr=descr.example,ratioEdges=0.5,
                 max.features=50, type=type,goParallel=goParallel,verbose=TRUE)
    
    
    
    #X<-rbind(trainD2C@X,testD2C@X)
    #Y<-c(trainD2C@Y,testD2C@Y)
    X<-trainD2C@X
    Y<-trainD2C@Y
    N=NROW(X)
    Itr=1:NROW(trainD2C@X) 
    Itr=sample(N,round(N/5))
    Its=setdiff(1:N,Itr)
    Xtr=X[Itr,2:50]
    Ytr=Y[Itr]
    Xts=X[Its,2:50]
    Yts=Y[Its]
    RF <- randomForest(x =Xtr ,y = factor(Ytr),importance=TRUE)
    Yhat=predict(RF,Xts)
    print(mean(Yhat!=Yts))
  }
  
  if (savefile)
    save(file=namefile,list=c("trainD2C","testDAG"))
}
##stopCluster(cl)

## number of DAGs used for testing
if (savefile)
  load(namefile)

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
  
  cat("Dim test dataset"=dim(observedData),"\n")
  
  ## inference of networks with bnlearn package
  
  ## Ahat.GS<-(amat(gs(data.frame(observedData))))
  Ahat.IAMB<-(amat(iamb(data.frame(observedData),alpha=0.01)))
  Ahat.PC<-(amat(si.hiton.pc(data.frame(observedData),alpha=0.01)))
  Ahat.GS<-Ahat.IAMB
  igraph.GS<-igraph::graph.adjacency(Ahat.GS)
  igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
  igraph.PC<-igraph::graph.adjacency(Ahat.PC)
  
  
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
    pred.D2C = predict(trainD2C,i,j, observedData)
    
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
  cat("\n r=",r," BER.D2C=",mean(BER.D2C), "BER.IAMB=",mean(BER.IAMB),
      "BER.GS=",mean(BER.GS),"BER.PC=",mean(BER.PC),"#0=",length(which(Ytrue==0))/length(Ytrue),"\n")
  
}

