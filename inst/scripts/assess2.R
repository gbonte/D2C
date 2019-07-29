rm(list=ls())
RNGversion("3.5.1")
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


file="traintestSTAR.1000.10.100.RData"
file2=paste("../D2Cdata2/",file,sep="")
resultsfile<-paste("../D2Cdata2/","ass3.",file,sep="")

load(file2)

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

print(testDAG@NDAG)
r=1
ncores=2
nedges=20
cl <- makeForkCluster(ncores)
registerDoParallel(cl)


while ( r <=testDAG@NDAG){
  
  FF<-foreach (i=r:(r+ncores-1)) %do%{
    ## for (i in  r:(r+ncores-1) ){
    set.seed(i)
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
    Yhat.GRA<-NULL
    Ytrue<-NULL
    
    observedData<-testDAG@list.observationsDAGs[[i]]
    
    N<-NROW(observedData)
    n<-NCOL(observedData)
    
    
    if (!any(is.nan(observedData))){
      trueDAG<-testDAG@list.DAGs[[i]]
      
      cat("Dim test dataset",dim(observedData),"\n")
      
      ## inference of networks with bnlearn package
      
      
      Ahat.IAMB<-(amat(iamb(data.frame(observedData[1:10,]),alpha=0.01)))
      
      print("IAMB")
      Ahat.IIAMB<-Ahat.IAMB
      
      print("IIAMB")
      Ahat.FIAMB<-Ahat.IAMB
      print("FIAMB")
      Ahat.PC<-Ahat.IAMB
      print("PC")
      Ahat.TABU<-Ahat.IAMB
      print("TABU")
      Ahat.RS<-Ahat.IAMB
      print("RSMAX")
      
      Ahat.GS<-Ahat.PC
      Ahat.PCalg<-Ahat.PC
      Ahat.KPCalg<-Ahat.PC
      
      
      
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
      set.seed(i)
      ## selection of a balanced subset of edges for the assessment
      Nodes=graph::nodes(trueDAG)
      max.edges<-min(nedges,length(gRbase::edgeList(trueDAG)))
      subset.edges = matrix(unlist(sample(gRbase::edgeList(trueDAG),size = max.edges,replace = F)),
                            ncol=2,byrow = TRUE)
      subset.edges = rbind(subset.edges,t(replicate(n =max.edges ,
                                                    sample(Nodes,size=2,replace = FALSE))))
      
      
      for(jj in 1:NROW(subset.edges)){
        i =as(subset.edges[jj,1],"numeric")
        j =as(subset.edges[jj,2],"numeric")
        
        
        
        Yhat.D2C.1<-c(Yhat.D2C.1,1)
        Yhat.D2C.2<-Yhat.D2C.1
        Yhat.D2C.3<-Yhat.D2C.1 
        prob.D2C.1<-c(prob.D2C.1,0)
        pred.GRANGER=1-grangertest(observedData[,i],observedData[,j], order = 2)$"Pr(>F)"[2]
        Yhat.GRA<-c(Yhat.GRA,as.numeric(pred.GRANGER>0.5) )
        
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
        
        
      } ## for jj
    }
    list(BER.D2C=BER(Ytrue,Yhat.D2C.1), AUC.D2C=AUC(Ytrue,prob.D2C.1),
         BER.IAMB=BER(Ytrue,Yhat.IAMB), BER.IIAMB=BER(Ytrue,Yhat.IIAMB),
         BER.FIAMB=BER(Ytrue,Yhat.FIAMB), 
         BER.GS=BER(Ytrue,Yhat.GS),
         BER.PC=BER(Ytrue,Yhat.PC),
         BER.TABU=BER(Ytrue,Yhat.TABU),
         BER.RS=BER(Ytrue,Yhat.RS),
         BER.PCalg=BER(Ytrue,Yhat.PCalg),
         BER.KPCalg=BER(Ytrue,Yhat.KPCalg),
         BER.GRA=BER(Ytrue,Yhat.GRA),
         N=N,n=n)
  } ## foreach
  
  for (i in 1:length(FF)){
    
    ## computation of Balanced Error Rate
    BER.D2C.1<-c(BER.D2C.1,FF[[i]]$BER.D2C)
    BER.D2C.2<-BER.D2C.1
    BER.D2C.3<-BER.D2C.2
    AUC.D2C.1<-c(AUC.D2C.1,FF[[i]]$AUC.D2C)
    BER.IAMB<-c(BER.IAMB,FF[[i]]$BER.IAMB)
    BER.IIAMB<-c(BER.IIAMB,FF[[i]]$BER.IIAMB)
    BER.FIAMB<-c(BER.FIAMB,FF[[i]]$BER.FIAMB)
    BER.GS<-c(BER.GS,FF[[i]]$BER.GS)
    BER.PC<-c(BER.PC,FF[[i]]$BER.PC)
    BER.TABU<-c(BER.TABU,FF[[i]]$BER.TABU)
    BER.RS<-c(BER.RS,FF[[i]]$BER.RS)
    BER.PCalg<-c(BER.PCalg,FF[[i]]$BER.PCalg)
    BER.KPCalg<-c(BER.KPCalg,FF[[i]]$BER.KPCalg)
    BER.GRA<-c(BER.GRA,FF[[i]]$BER.GRA)
    ALLN<-c(ALLN,FF[[i]]$N)
    ALLn<-c(ALLn,FF[[i]]$n)
  }
  
  cat("\n r=",r," AUC.D2C.1=",mean(AUC.D2C.1)," BER.D2C.1=",mean(BER.D2C.1),
      " BER.D2C.2=",mean(BER.D2C.2)," BER.D2C.3=",mean(BER.D2C.3),
      "BER.IAMB=",mean(BER.IAMB),"BER.IIAMB=",mean(BER.IIAMB),"BER.FIAMB=",mean(BER.FIAMB),
      "BER.GS=",mean(BER.GS),"BER.PC=",mean(BER.PC),
      "BER.PCalg=",mean(BER.PCalg),"BER.TABU=",mean(BER.TABU),
      "BER.RS=",mean(BER.RS),
      "BER.KPCalg=",mean(BER.KPCalg),
      "BER.GRA=",mean(BER.GRA),"\n")
  r=r+ncores
  #save(file=resultsfile,list=c("BER.D2C.1","BER.IAMB","BER.FIAMB","BER.IIAMB",
  #                             "BER.GS","BER.PC","BER.TABU","BER.RS","AUC.D2C.1",
  #                             "BER.PCalg","BER.KPCalg","BER.GRA",
  #                             "ALLN","ALLn"))
} ## while







