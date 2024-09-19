### D2C Code related to the D2C algorithm 
## from "From Dependency to Causality: A Machine Learning Approach" in http://jmlr.org/papers/v16/bontempi15a.html
## 
## Gianluca Bontempi, mlg.ulb.ac.be

## Small scale assessment of D2C wrt state of the art "bnlearn" inference algorithms 

rm(list=ls())
library(D2C)
require(RBGL)
require(gRbase)
require(igraph)
require(graph)
require(pcalg)
 

set.seed(0)


###########################


noNodes<-c(10,20)
## range of number of nodes

N<-c(50,100)
## range of number of samples

sd.noise<-c(0.1,0.25)
## range of values for standard deviation of additive noise 

NDAG=10
## number of DAGs to be created and simulated

cat(paste("\n Generating", NDAG, " training DAGs ... "))
type="is.parent"
trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic"), 
              seed=1,sdn=sd.noise,additive=c(TRUE,FALSE),verbose=FALSE,
              maxV=3,weights=c(0.5,1))



###########################
cat("\n Computing descriptors ... ")
descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=5,acc=TRUE,lin=FALSE,residual=TRUE)

D2C<-new("D2C",sDAG=trainDAG, npar=10,
         descr=descr.example,ratioEdges=0.1,max.features=20,type=type,
         verbose=FALSE)

###########################
cat("\n Learning classifier ... ")
trainD2C<-makeModel(D2C,classifier="RF",EErep=15,verbose=FALSE)

###########################
cat("\n Generating testing DAGs ... ")
NDAG.test=10


testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
             functionType = c("linear","quadratic","sigmoid"), 
             quantize=c(TRUE,FALSE),
             seed=101,sdn=sd.noise,verbose=FALSE)


###########################
cat("\n D2C assessment in terms of BER (Balanced Error Rate)... ")

if (!require(bnlearn)){
  install.packages("bnlearn", repos="http://cran.rstudio.com/")
  library(bnlearn)
}
gendata=FALSE
Yhat.D2C<-NULL
Yhat.IAMB<-NULL
Yhat.GS<-NULL
Ytrue<-NULL
for (r in 1:testDAG@NDAG){
  set.seed(r)
  if (gendata){
    g<-gendataDAG(100,sample(noNodes[1]:noNodes[2],1))
    trueDAG<-g$DAG
    observedData <- g$data
  } else {
    observedData<-testDAG@list.observationsDAGs[[r]]
    trueDAG<-testDAG@list.DAGs[[r]]
  }
  if (NROW(observedData)>20){
    #
    
    ## inference of networks with bnlearn package
    Ahat.GS<-amat(gs(data.frame(observedData),alpha=0.01))
    Ahat.IAMB<-(amat(iamb(data.frame(observedData),alpha=0.01)))
    
    graphTRUE<- as.adjMAT(trueDAG)
    igraphTRUE<-graph.adjacency(as(trueDAG,"matrix"))
    
    colnames(Ahat.IAMB)=colnames(graphTRUE)
    rownames(Ahat.IAMB)=rownames(graphTRUE)
    igraphIAMB<-graph.adjacency(Ahat.IAMB)
    
    colnames(Ahat.GS)=colnames(graphTRUE)
    rownames(Ahat.GS)=rownames(graphTRUE)
    igraphGS<-graph.adjacency(Ahat.GS)
    
    
    ## selection of a balanced subset of edges for the assessment
    Nodes=graph::nodes(trueDAG)
    max.edges<-min(10,length(edgeList(as(trueDAG,'matrix'))))
    subset.edges = matrix(unlist(sample(edgeList(as(trueDAG,'matrix')),
                                        size = max.edges,replace = F)),ncol=2,byrow = TRUE)
    subset.edges = rbind(subset.edges,t(replicate(n =max.edges,
                                                  sample(Nodes,size=2,replace = FALSE))))
    
    
    for(jj in 1:NROW(subset.edges)){
      i =as(subset.edges[jj,1],"numeric");
      j =as(subset.edges[jj,2],"numeric") ;
      pred.D2C = predict(trainD2C,i,j, observedData,rep=3)
      
      Yhat.D2C<-c(Yhat.D2C,as.numeric(pred.D2C$response))
     
      
      Yhat.IAMB<-c(Yhat.IAMB,is.what(igraphIAMB,
                                     subset.edges[jj,1],subset.edges[jj,2],type))
      
      Yhat.GS<-c(Yhat.GS,is.what(igraphGS,
                                 subset.edges[jj,1],subset.edges[jj,2],type))
      Ytrue<-c(Ytrue,is.what(igraphTRUE,
                             subset.edges[jj,1],subset.edges[jj,2],type))
      ## computation of Balanced Error Rate
      BER.D2C<-D2C::BER(Ytrue,Yhat.D2C)
      BER.GS<-D2C::BER(Ytrue,Yhat.GS)
      BER.IAMB<-D2C::BER(Ytrue,Yhat.IAMB)
      
      ## ----echo=TRUE, eval=FALSE-----------------------------------------------
      cat("\n  Test DAG", r, ",", length(Ytrue), 
          "edges tested: BER.D2C=",BER.D2C, "BER.IAMB=",BER.IAMB,"BER.GS=",BER.GS,"\n")
    }
    
  }
  
}
