rm(list=ls())
library(D2C)
require(bnlearn)
library(pcalg)
library(kpcalg)
library(lmtest)
library(doParallel)
type="is.parent"



set.seed(0)


goParallel=TRUE
savefile<-TRUE
namefile<-"./data/traintestTSERIES.RData"
load(namefile)
noNodes<-6 ## max lag 

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
MXSX=3
print(testDAG@NDAG)

for ( r in 1:testDAG@NDAG){
  set.seed(r)
  observedData<-testDAG@list.observationsDAGs[[r]]
  n<-NCOL(observedData)
  trueDAG<-testDAG@list.DAGs[[r]]
  
  
  cat("Dim test dataset"=dim(observedData),"\n")
  
  ## inference of networks with bnlearn package
  Ahat.IAMB<-(amat(iamb(data.frame(observedData),
                        alpha=0.01,max.sx=MXSX)))
  print("Done IAMB")
  Ahat.PC<-(amat(si.hiton.pc(data.frame(observedData),alpha=0.01,max.sx=MXSX)))
  print("Done Hiton PC")
  
  if (n<1000){
    Ahat.GS<-(amat(gs(data.frame(observedData),max.sx=MXSX)))
    
    suffStat <- list(C = cor(observedData),n=NROW(observedData))
    normal.pag <- pc(suffStat, indepTest=gaussCItest, alpha = 0.01, ,m.max=MXSX,
                     verbose=FALSE,p=NCOL(observedData),
                     numCores=3)
    Ahat.GS<-as(normal.pag, "matrix") 
    print("Done pc")

    igraph.GS<-igraph::graph.adjacency(Ahat.GS)
    igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
    igraph.PC<-igraph::graph.adjacency(Ahat.PC)

    
    graphTRUE<- gRbase::as.adjMAT(trueDAG)
    igraph.TRUE<-igraph::graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),as.character(1:NCOL(graphTRUE))])
    
    ## selection of a balanced subset of edges for the assessment
    Nodes=graph::nodes(trueDAG)
    max.edges<-min(25,length(gRbase::edgeList(trueDAG)))
    subset.edges = matrix(unlist(sample(gRbase::edgeList(trueDAG),size = max.edges,replace = F)),
                          ncol=2,byrow = TRUE)
    subset.edges = rbind(subset.edges,t(replicate(n =max.edges ,
                                                  sample(Nodes,size=2,replace = FALSE))))
    
    
    for(jj in 1:NROW(subset.edges)){
      i =as(subset.edges[jj,1],"numeric"); # parent
      j =as(subset.edges[jj,2],"numeric") ; ## child
      
      fs<-timecauses(NCOL(observedData),noNodes,j)
      if (TRUE){
        
        pred.D2C = predict(trainD2C,i,j, observedData,rep=3)
    
        Yhat.D2C<-c(Yhat.D2C,as.numeric(pred.D2C$response))
        
        Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,i,j,"is.parent"))
        Yhat.GS<-c(Yhat.GS,is.what(igraph.GS,i,j,"is.parent"))
        Yhat.PC<-c(Yhat.PC,is.what(igraph.PC,i,j,"is.parent"))
       
        Ytrue<-c(Ytrue,is.what(igraph.TRUE,i,j,"is.parent")) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
        
        cat(".")
      }
    }
    
    
    ## computation of Balanced Error Rate
    BER.D2C<-BER(Ytrue,Yhat.D2C)
    
    BER.IAMB<-BER(Ytrue,Yhat.IAMB)
    BER.GS<-BER(Ytrue,Yhat.GS)
    BER.PC<-BER(Ytrue,Yhat.PC)
    
    
    cat("\n r=",r," BER.D2C=",mean(BER.D2C), 
        "BER.IAMB=",mean(BER.IAMB),
        "BER.GS=",mean(BER.GS),"BER.PC=",mean(BER.PC),
        "#0=",length(which(Ytrue==0))/length(Ytrue),"\n")
    
  }
}

