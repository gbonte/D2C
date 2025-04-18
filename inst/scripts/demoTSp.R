rm(list=ls())
library(D2C)
require(bnlearn)
library(pcalg)
#library(kpcalg)
library(lmtest)
library(doParallel)
type="is.ancestor"

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
  if (type=="is.distance")
    return(dagdistance(iDAG,i,j))
}

noNodes<-c(7,8)
## range of number of nodes

N<-c(150,200)
## range of number of samples

NDAG=50
## number of DAGs to be created and simulated
NDAG.test=500
nseries=5
sdev<-c(0.1,0.2)

goParallel=TRUE
savefile<-TRUE
namefile<-"../data/demoTS2.RData"
if (goParallel){
  ncores=10
  cl <- makeForkCluster(ncores)
  registerDoParallel(cl)
  dN=3*ncores
}




trainDAG<-new("simulatedTS",NDAG=NDAG, N=N, noNodes=noNodes,
              seed=10,sdn=sdev,goParallel=goParallel,
              nseries=nseries,typeser=c(15:23))


descr.example<-new("D2C.descriptor",bivariate=TRUE,
                   ns=5,maxs=10,acc=TRUE,lin=TRUE,
                   struct=TRUE, boot="mimr")


trainD2C<-new("D2C",sDAG=trainDAG,
              descr=descr.example,ratioEdges=1,
              max.features=20, type=type,goParallel=goParallel,
              verbose=TRUE)

print(NROW(trainD2C@origX))

trainD2C.1<-makeModel(trainD2C,classifier="RF",EErep=2)

testDAG<-new("simulatedTS",NDAG=NDAG.test, N=N, noNodes=noNodes,
             seed=101,sdn=sdev,goParallel=goParallel,
             nseries=nseries,
             typeser=1:14)

if (savefile)
  save(file=namefile,list=c("trainD2C.1","testDAG"))

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
    igraph.GS<-igraph::graph.adjacency(Ahat.GS)
    igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
    igraph.PC<-igraph::graph.adjacency(Ahat.PC)
    igraph.KPC<-igraph::graph.adjacency(Ahat.KPC)
    
    graphTRUE<- gRbase::as.adjMAT(trueDAG)
    igraph.TRUE<-igraph::graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),as.character(1:NCOL(graphTRUE))])
    
    ## selection of a balanced subset of edges for the assessment
    Nodes=graph::nodes(trueDAG)
    max.edges<-min(100,length(gRbase::edgeList(trueDAG)))
    subset.edges = matrix(unlist(sample(gRbase::edgeList(trueDAG),
                                        size = max.edges,replace = F)),
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
      Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,i,j,type))
      Yhat.GS<-c(Yhat.GS,is.what(igraph.GS,i,j,type))
      Yhat.PC<-c(Yhat.PC,is.what(igraph.PC,i,j,type))
      Yhat.KPC<-c(Yhat.KPC,is.what(igraph.KPC,i,j,type))
      
      Ytrue<-c(Ytrue,is.what(igraph.TRUE,i,j,type)) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
      
      cat(".")
      
    }
    if (type=="is.distance"){
      
      BER.D2C<-mean((Ytrue-Yhat.D2C)^2)/var(Ytrue)
      BER.D2C.1<-mean((Ytrue-Yhat.D2C.1)^2)/var(Ytrue)
      BER.IAMB<-mean((Ytrue-Yhat.IAMB)^2)/var(Ytrue)
      BER.GS<-mean((Ytrue-Yhat.GS)^2)/var(Ytrue)
      BER.PC<-mean((Ytrue-Yhat.PC)^2)/var(Ytrue)
      BER.KPC<-mean((Ytrue-Yhat.KPC)^2)/var(Ytrue)
      BER.GRA<-mean((Ytrue-Yhat.GRA)^2)/var(Ytrue)
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
    }
    cat("\n r=",r," BER.D2C=",mean(BER.D2C), " BER.D2C.1=",mean(BER.D2C.1),
        "BER.IAMB=",mean(BER.IAMB),"BER.GRA=",mean(BER.GRA),
        "BER.GS=",mean(BER.GS),"BER.PC=",mean(BER.PC),"BER.KPC=",mean(BER.KPC),
        "#0=",length(which(Ytrue==0))/length(Ytrue),"\n")
    
  }
}


