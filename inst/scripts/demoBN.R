## ----echo=TRUE,results=TRUE----------------------------------------------
rm(list=ls())
library(D2C)
require(RBGL)
require(gRbase)
require(igraph)
require(graph)
require(pcalg)
 

set.seed(0)

noNodes<-c(5,20)
## range of number of nodes

N<-c(50,100)
## range of number of samples

sd.noise<-c(0.1,0.25)
## range of values for standard deviation of additive noise 

NDAG=50
## number of DAGs to be created and simulated

type="is.parent"
trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid"), 
              seed=1,sdn=0,additive=c(TRUE,FALSE),verbose=TRUE,maxV=3,weights=c(0.5,1))



## ---- echo=TRUE----------------------------------------------------------
print(trainDAG@NDAG)

## ---- echo=TRUE----------------------------------------------------------
print(trainDAG@list.DAGs[[1]])

print(dim(trainDAG@list.observationsDAGs[[1]]))

## ----echo=TRUE,eval=TRUE-------------------------------------------------
descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=5,acc=TRUE,lin=FALSE,residual=TRUE)

D2C<-new("D2C",sDAG=trainDAG, npar=10,
         descr=descr.example,ratioEdges=0.15,max.features=30,type=type,
         verbose=TRUE)


trainD2C<-makeModel(D2C,classifier="RF",EErep=15)

## ----echo=TRUE,results=FALSE---------------------------------------------
print(dim(trainD2C@X))
print(table(trainD2C@Y))

## ----echo=TRUE,results=FALSE---------------------------------------------
print(trainD2C@mod)


## ----echo=TRUE-----------------------------------------------------------
NDAG.test=10


testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
             functionType = c("linear","quadratic","sigmoid"), 
             quantize=c(TRUE,FALSE),
             seed=101,sdn=sd.noise,verbose=TRUE)

## ----echo=TRUE,results=FALSE,message=FALSE, warning=FALSE,eval=FALSE-----

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
    max.edges<-min(10,length(edgeList(trueDAG)))
    subset.edges = matrix(unlist(sample(edgeList(trueDAG),
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
      cat("\n  DAG", r, ",", length(Ytrue), 
          " edges tested: BER.D2C=",BER.D2C, "BER.IAMB=",BER.IAMB,"BER.GS=",BER.GS,"\n")
    }
    
  }
  
}
