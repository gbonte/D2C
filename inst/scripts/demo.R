## ----echo=TRUE,results=TRUE----------------------------------------------
rm(list=ls())
library(D2C)
require(RBGL)
require(gRbase)
require(igraph)
require(graph)

set.seed(0)

noNodes<-c(100,150)
## range of number of nodes

N<-c(100,200)
## range of number of samples

sd.noise<-c(0.2,0.5)
## range of values for standard deviation of additive noise 

NDAG=50
## number of DAGs to be created and simulated


trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid"), 
              seed=0,sdn=sd.noise,quantize=c(TRUE,FALSE),verbose=TRUE)



## ---- echo=TRUE----------------------------------------------------------
print(trainDAG@NDAG)

## ---- echo=TRUE----------------------------------------------------------
print(trainDAG@list.DAGs[[1]])
print(dim(trainDAG@list.observationsDAGs[[1]]))

## ----echo=TRUE,eval=TRUE-------------------------------------------------
descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=5,acc=TRUE,lin=FALSE)

D2C<-new("D2C",sDAG=trainDAG, 
              descr=descr.example,ratioEdges=0.1,max.features=30,verbose=TRUE)


trainD2C<-makeModel(D2C,classifier="RF")

## ----echo=TRUE,results=FALSE---------------------------------------------
print(dim(trainD2C@X))
print(table(trainD2C@Y))

## ----echo=TRUE,results=FALSE---------------------------------------------
print(trainD2C@mod)


## ----echo=TRUE-----------------------------------------------------------
NDAG.test=100


testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
             functionType = c("linear","quadratic","sigmoid"), 
             quantize=c(TRUE,FALSE),
             seed=101,sdn=sd.noise,verbose=TRUE)

## ----echo=TRUE,results=FALSE,message=FALSE, warning=FALSE,eval=FALSE-----
 
  if (!require(bnlearn)){
    install.packages("bnlearn", repos="http://cran.rstudio.com/")
    library(bnlearn)
    }
Yhat.D2C<-NULL
Yhat.IAMB<-NULL
Yhat.GS<-NULL
Ytrue<-NULL
  for (r in 1:testDAG@NDAG){
    set.seed(r)
    observedData<-testDAG@list.observationsDAGs[[r]]
    trueDAG<-testDAG@list.DAGs[[r]]
  
    ## inference of networks with bnlearn package
    Ahat.GS<-amat(gs(data.frame(observedData),alpha=0.01))
    Ahat.IAMB<-(amat(iamb(data.frame(observedData),alpha=0.01)))
  
    graphTRUE<- as.adjMAT(trueDAG)
  
  
    ## selection of a balanced subset of edges for the assessment
    Nodes=nodes(trueDAG)
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
      Yhat.IAMB<-c(Yhat.IAMB,Ahat.IAMB[i,j])
      Yhat.GS<-c(Yhat.GS,Ahat.GS[i,j])
      Ytrue<-c(Ytrue,graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
      ## computation of Balanced Error Rate
      BER.D2C<-BER(Ytrue,Yhat.D2C)
      BER.GS<-BER(Ytrue,Yhat.GS)
      BER.IAMB<-BER(Ytrue,Yhat.IAMB)
      
      ## ----echo=TRUE, eval=FALSE-----------------------------------------------
      cat("\n  DAG", r, ",", length(Ytrue), " edges tested: BER.D2C=",BER.D2C, "BER.IAMB=",BER.IAMB,"BER.GS=",BER.GS,"\n")
      }
  
  
 
}
## ----echo=TRUE,eval=FALSE------------------------------------------------
#  data(alarm)
#  
#  graphTRUE<-true.net
#  set.seed(0)
#  
#  observedData<-dataset
#  
#  w.const<-which(apply(observedData,2,sd)<0.1)
#  if (length(w.const)>0){
#    observedData<-observedData[,-w.const]
#    graphTRUE<-graphTRUE[-w.const,-w.const]
#    }
#  indn<-sort(apply(graphTRUE,1,sum)+apply(graphTRUE,2,sum),decr=TRUE,index=TRUE)$ix[1:100]
#  observedData<-observedData[,indn]
#  graphTRUE<-graphTRUE[indn,indn]
#  

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  n<-NCOL(observedData)
#  
#  
#  Ahat.GS<-amat(gs(data.frame(observedData)))
#  Ahat.IAMB<-amat(iamb(data.frame(observedData),alpha=0.05))

## ----echo=TRUE,results=FALSE,message=FALSE, warning=FALSE,eval=FALSE-----
#  Yhat.D2C<-NULL
#  Yhat.GS<-NULL
#  Yhat.IAMB<-NULL
#  Ytrue<-NULL
#  
#  
#  for (i in 1:n){
#    ## creation of a balanced test set
#    ind1<-which(graphTRUE[i,]==1)
#    ind0<-setdiff(setdiff(1:n,i),ind1)
#    ind<-c(ind1,ind0[1:length(ind1)])
#  
#    FF<-foreach(j=ind) %do%{
#       list(Yhat=as.numeric(predict(trainD2C,i,j,observedData)$response)  -1,
#           Yhat2=Ahat.GS[i,j],Yhat3=Ahat.IAMB[i,j],Ytrue=graphTRUE[i,j])
#      }
#    Yhat.D2C<-c(Yhat.D2C,unlist(lapply(FF,"[[",1) ))
#    Yhat.GS<-c(Yhat.GS,unlist(lapply(FF,"[[",2)))
#    Yhat.IAMB<-c(Yhat.IAMB,unlist(lapply(FF,"[[",3)))
#    Ytrue<-c(Ytrue,unlist(lapply(FF,"[[",4)))
#    }

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  cat("\n BER.D2C=",BER(Ytrue,Yhat.D2C), "BER.GS=",BER(Ytrue,Yhat.GS),
#      "BER.IAMB=",BER(Ytrue,Yhat.IAMB),
#      "\n")

