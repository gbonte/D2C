rm(list=ls())
#RNGversion("3.5.1")
library(D2C)
require(bnlearn)
library(pcalg)
library(kpcalg)
library(ROCR)
library(foreach)
library(doParallel)
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


## range of number of nodes


file="trainD2Cp.5000.50.is.parent.RData"
filetrain=paste("./data/",file,sep="")
resultsfile<-paste("./data/","ass.",file,sep="")

load(filetrain)

file="testDAG.100.100.RData"
filetest=paste("./data/",file,sep="")
load(filetest)

BER.D2C.1<-NULL
AUC.D2C.1<-NULL
BER.D2C.2<-NULL
BER.D2C.3<-NULL
BER.IAMB<-NULL
BER.IIAMB<-NULL
BER.FIAMB<-NULL
BER.GS<-NULL
BER.PC<-NULL
BER.TABU<-NULL
BER.RS<-NULL
BER.PCalg<-NULL
BER.KPCalg<-NULL
BER.GRA<-NULL
ALLN<-NULL
ALLn<-NULL
cat("NROW(X)=",NROW(trainD2C@origX),"\n")
print(testDAG@NDAG)
r=1
ncores=2
nedges=40
cl <- makeForkCluster(ncores)
registerDoParallel(cl)



##FF<-foreach (i=r:(r+ncores-1)) %do%{
prob.D2C.1<-NULL
Yhat.D2C.1<-NULL
Yhat.D2C.2<-NULL
Yhat.D2C.3<-NULL
Yhat.IAMB<-NULL
Yhat.IIAMB<-NULL
Yhat.FIAMB<-NULL
Yhat.GS<-NULL
Yhat.PC<-NULL
Yhat.TABU<-NULL
Yhat.RS<-NULL
Yhat.PCalg<-NULL
Yhat.KPCalg<-NULL
Ytrue<-NULL

for (ii in  1:length(testDAG@list.observationsDAGs) ){
  set.seed(ii)
  observedData<-testDAG@list.observationsDAGs[[ii]]
  
  N<-NROW(observedData)
  n<-NCOL(observedData)
  
  MX.SX=4
  
  if (!any(is.nan(observedData))){
    trueDAG<-testDAG@list.DAGs[[ii]]
    
    
    if (remNodes>0){
      n=NCOL(observedData)
      torem=sample(1:n,remNodes)
      toremch=as.character(torem)
      for (r in 1:remNodes)
        trueDAG<-graph::removeNode(toremch[r],trueDAG)
      observedData=observedData[,-torem]
      graph::nodes(trueDAG)<-as.character(1:NCOL(observedData))
      
    }
    
    cat("Dim test dataset",dim(observedData),"\n")
    
    ## inference of networks with bnlearn package
    
    
    Ahat.IAMB<-(amat(iamb(data.frame(observedData),alpha=0.01,max.sx=MX.SX)))
    
    print("IAMB")
    Ahat.IIAMB<-(amat(inter.iamb(data.frame(observedData),alpha=0.01,,max.sx=MX.SX)))
    print("IIAMB")
    Ahat.FIAMB<-(amat(fast.iamb(data.frame(observedData),alpha=0.01,,max.sx=MX.SX)))
    print("FIAMB")
    Ahat.PC<-(amat(si.hiton.pc(data.frame(observedData),alpha=0.01,,max.sx=MX.SX)))
    print("PC")
    Ahat.TABU<-(amat(tabu(data.frame(observedData),max.iter = 1000)))
    print("TABU")
    Ahat.RS<-Ahat.TABU 
    print("RSMAX")
    if (n<20){
      Ahat.RS<-(amat(rsmax2(data.frame(observedData))))
      Ahat.GS<-(amat(gs(data.frame(observedData))))
      print("GS")
      suffStat <- list(C = cor(observedData),n=NROW(observedData))
      normal.pag <- pc(suffStat, indepTest=gaussCItest, alpha = 0.01, ,m.max=MX.SX,
                       verbose=FALSE,p=NCOL(observedData),
                       numCores=3)
      Ahat.PCalg<-as(normal.pag, "matrix") 
      print("Done pcalg")
      # kpag <- kpc(suffStat = list(data=observedData[sample(N,min(N,100)),], 
      #                            ic.method="hsic.gamma"),
      #            indepTest=kernelCItest, alpha = 0.01, skel.method = "stable.fast",
      #            verbose=FALSE,p=NCOL(observedData),m.max=2)
      Ahat.KPCalg<-Ahat.PCalg 
      #browser()
      print("Done kpc")
    } else {
      Ahat.GS<-Ahat.PC
      Ahat.PCalg<-Ahat.PC
      Ahat.KPCalg<-Ahat.PC
      
    }
    
    igraph.PCalg<-igraph::graph.adjacency(Ahat.PCalg)
    igraph.KPCalg<-igraph::graph.adjacency(Ahat.KPCalg)
    igraph.GS<-igraph::graph.adjacency(Ahat.GS)
    igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
    igraph.IIAMB<-igraph::graph.adjacency(Ahat.IIAMB)
    igraph.FIAMB<-igraph::graph.adjacency(Ahat.FIAMB)
    igraph.PC<-igraph::graph.adjacency(Ahat.PC)
    igraph.TABU<-igraph::graph.adjacency(Ahat.TABU)
    igraph.RS<-igraph::graph.adjacency(Ahat.RS)
    
    
    graphTRUE<- gRbase::as.adjMAT(trueDAG)
    igraph.TRUE<-igraph::graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),
                                                   as.character(1:NCOL(graphTRUE))])
    
    ## selection of a balanced subset of edges for the assessment
    Nodes=graph::nodes(trueDAG)
    max.edges<-min(nedges,length(gRbase::edgeList(trueDAG)))
    subset.edges = matrix(unlist(sample(gRbase::edgeList(trueDAG),
                                        size = max.edges,replace = F)),
                          ncol=2,byrow = TRUE)
    subset.edges = rbind(subset.edges,t(replicate(n =max.edges ,
                                                  sample(Nodes,size=2,replace = FALSE))))
    
    for(jj in 1:NROW(subset.edges)){
      i =as(subset.edges[jj,1],"numeric")
      j =as(subset.edges[jj,2],"numeric")
  
      pred.D2C = predict(trainD2C,i,j, observedData,rep=5)
      Yhat.D2C.1<-c(Yhat.D2C.1,pred.D2C$response)
      Yhat.D2C.2<-Yhat.D2C.1
      Yhat.D2C.3<-Yhat.D2C.1 
      prob.D2C.1<-c(prob.D2C.1,pred.D2C$prob)
      
      Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,i,j))
      Yhat.IIAMB<-c(Yhat.IIAMB,is.what(igraph.IIAMB,i,j))
      Yhat.FIAMB<-c(Yhat.FIAMB,is.what(igraph.FIAMB,i,j))
      Yhat.GS<-c(Yhat.GS,is.what(igraph.GS,i,j))
      Yhat.PC<-c(Yhat.PC,is.what(igraph.PC,i,j))
      Yhat.TABU<-c(Yhat.TABU,is.what(igraph.TABU,i,j))
      Yhat.RS<-c(Yhat.RS,is.what(igraph.RS,i,j))
      Yhat.PCalg<-c(Yhat.PCalg,is.what(igraph.PCalg,i,j))
      Yhat.KPCalg<-c(Yhat.KPCalg,is.what(igraph.KPCalg,i,j))
      Ytrue<-c(Ytrue,is.what(igraph.TRUE,i,j)) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
      
      cat(".")
    } ## for jj
    BER.D2C=D2C::BER(Ytrue,Yhat.D2C.1)
    AUC.D2C=D2C::AUC(Ytrue,prob.D2C.1)
    BER.IAMB=D2C::BER(Ytrue,Yhat.IAMB)
    BER.IIAMB=D2C::BER(Ytrue,Yhat.IIAMB)
    BER.FIAMB=D2C::BER(Ytrue,Yhat.FIAMB) 
    BER.GS=D2C::BER(Ytrue,Yhat.GS)
    BER.PC=D2C::BER(Ytrue,Yhat.PC)
    BER.TABU=D2C::BER(Ytrue,Yhat.TABU)
    BER.RS=D2C::BER(Ytrue,Yhat.RS)
    BER.PCalg=D2C::BER(Ytrue,Yhat.PCalg)
    BER.KPCalg=D2C::BER(Ytrue,Yhat.KPCalg)
    
    
    cat("\n ii=",ii," AUC.D2C=",mean(AUC.D2C)," BER.D2C=",mean(BER.D2C),
       "BER.IAMB=",mean(BER.IAMB),"BER.IIAMB=",mean(BER.IIAMB),"BER.FIAMB=",mean(BER.FIAMB),
        "BER.GS=",mean(BER.GS),"BER.PC=",mean(BER.PC),
        "BER.PCalg=",mean(BER.PCalg),"BER.TABU=",mean(BER.TABU),
        "BER.RS=",mean(BER.RS),
        "BER.KPCalg=",mean(BER.KPCalg),
        "\n")
    
    
  }
  
}  
  
  
  
  
  