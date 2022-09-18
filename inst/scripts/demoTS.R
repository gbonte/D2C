rm(list=ls())
library(D2C)
require(bnlearn)
library(pcalg)
#library(kpcalg)
library(lmtest)
type="is.ancestor"

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
  if (type=="is.distance")
    return(dagdistance(iDAG,i,j))
}

noNodes<-c(15,50)
## range of number of nodes

N<-c(50,100)
## range of number of samples

NDAG=50
## number of DAGs to be created and simulated
NDAG.test=100
nseries=3
sdev<-c(0.1,0.3)

goParallel=FALSE
savefile<-FALSE
namefile<-"../D2Cdata2/traintestSTAR.200.100.RData"
if (TRUE){
  
  
  
  trainDAG<-new("simulatedTS",NDAG=NDAG, N=N, noNodes=noNodes,
                seed=10,sdn=sdev,goParallel=goParallel,nseries=nseries,typeser=c(15:23))
  
  
  descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=5,maxs=10,acc=TRUE,lin=FALSE,
                     struct=FALSE, boot="rank")
  
  
  trainD2C<-new("D2C",sDAG=trainDAG,
                descr=descr.example,ratioEdges=0.5,
                max.features=20, type=type,goParallel=goParallel,
                
                verbose=TRUE)
  

  trainD2C.1<-trainD2C
  print(NROW(trainD2C.1@origX))
  print(NROW(trainD2C@origX))
 
  trainD2C.1<-makeModel(trainD2C.1,classifier="RF",EErep=2)
  
  testDAG<-new("simulatedTS",NDAG=NDAG.test, N=N, noNodes=noNodes,
               seed=101,sdn=sdev,goParallel=goParallel,nseries=round(nseries/2),
               typeser=1:13)
  
  
  
  if (savefile)
    save(file=namefile,list=c("trainD2C","testDAG"))
}
##stopCluster(cl)

## number of DAGs used for testing
if (savefile)
  load(namefile)


print(colnames(trainD2C.1@X))
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
  n<-NCOL(observedData)
  trueDAG<-testDAG@list.DAGs[[r]]
  
  cat("Dim test dataset"=dim(observedData),"\n")
  
  ## inference of networks with bnlearn package
  Ahat.IAMB<-(amat(iamb(data.frame(observedData),#optimized=TRUE,
                        alpha=0.01,max.sx=3)))
  print("Done IAMB")
  Ahat.PC<-(amat(si.hiton.pc(data.frame(observedData),alpha=0.01)))
  
  if (n<100){
    Ahat.GS<-(amat(gs(data.frame(observedData))))
    
    
    
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
    graphTRUE<- gRbase::as.adjMAT(trueDAG)
    igraph.TRUE<-igraph::graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),as.character(1:NCOL(graphTRUE))])
    
    colnames(Ahat.IAMB)=colnames(graphTRUE)
    rownames(Ahat.IAMB)=rownames(graphTRUE)
    colnames(Ahat.GS)=colnames(graphTRUE)
    rownames(Ahat.GS)=rownames(graphTRUE)
    colnames(Ahat.PC)=colnames(graphTRUE)
    rownames(Ahat.PC)=rownames(graphTRUE)
    colnames(Ahat.KPC)=colnames(graphTRUE)
    rownames(Ahat.KPC)=rownames(graphTRUE)
    igraph.GS<-igraph::graph.adjacency(Ahat.GS)
    igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
    igraph.PC<-igraph::graph.adjacency(Ahat.PC)
    igraph.KPC<-igraph::graph.adjacency(Ahat.KPC)
    
    
    ## selection of a balanced subset of edges for the assessment
    Nodes=graph::nodes(trueDAG)
    max.edges<-min(100,length(gRbase::edgeList(trueDAG)))
    subset.edges = matrix(unlist(sample(gRbase::edgeList(trueDAG),size = max.edges,replace = F)),
                          ncol=2,byrow = TRUE)
    subset.edges = rbind(subset.edges,t(replicate(n =max.edges ,
                                                  sample(Nodes,size=2,replace = FALSE))))
    
    
    for(jj in 1:NROW(subset.edges)){
      i =as(subset.edges[jj,1],"numeric");
      j =as(subset.edges[jj,2],"numeric") ;
      pred.D2C = predict(trainD2C.1,i,j, observedData)
      pred.D2C.1 = pred.D2C #predict(trainD2C.1,i,j, observedData)
      pred.GRANGER=1-grangertest(observedData[,i],observedData[,j], order = 1)$"Pr(>F)"[2]
      Yhat.GRA<-c(Yhat.GRA,as.numeric(pred.GRANGER>0.5) )
      Yhat.D2C<-c(Yhat.D2C,pred.D2C$response)
      Yhat.D2C.1<-c(Yhat.D2C.1,pred.D2C.1$response)
      Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,i,j))
      Yhat.GS<-c(Yhat.GS,is.what(igraph.GS,i,j))
      Yhat.PC<-c(Yhat.PC,is.what(igraph.PC,i,j))
      Yhat.KPC<-c(Yhat.KPC,is.what(igraph.KPC,i,j))
      Ytrue<-c(Ytrue,is.what(igraph.TRUE,i,j)) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
      
      cat(".")
    }
    
    if (type=="is.distance"){
      browser()
    }
    else {
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
  }
}

