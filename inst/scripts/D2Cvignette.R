
rm(list=ls())
library(D2C)
require(foreach)
require(doParallel)
require(RBGL)
require(gRbase)
require(igraph)
require(graph)

set.seed(0)

noNodes<-c(50,100)
## range of number of nodes

N<-c(100,200)
##  range of number of samples

sd.noise<-c(0.2,0.5)
##  range of values for standard deviation of additive noise

NDAG=20
##  number of DAGs to be created and simulated


trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid"),
              seed=0,sdn=sd.noise,quantize=c(TRUE,FALSE),verbose=TRUE)



print(trainDAG@NDAG) 


print(trainDAG@list.DAGs[[1]])
print(dim(trainDAG@list.observationsDAGs[[1]]))


descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=5,acc=TRUE,lin=FALSE)

trainD2C<-new("D2C",sDAG=trainDAG,
              descr=descr.example,ratioEdges=0.1,max.features=30,verbose=TRUE)


trainD2C<-makeModel(trainD2C,classifier="RF",EErep=2)



print(dim(trainD2C@origX))
print(table(trainD2C@Y))


NDAG.test=10


testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
             functionType = c("linear","quadratic","sigmoid"),
             quantize=c(TRUE,FALSE),
             seed=101,sdn=sd.noise,verbose=TRUE)



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
   
   ##  inference of networks with bnlearn package
   Ahat.GS<-amat(gs(data.frame(observedData),alpha=0.01,max.sx=3))
   Ahat.IAMB<-(amat(iamb(data.frame(observedData),alpha=0.01,max.sx=3)))
   
   graphTRUE<- as.adjMAT(trueDAG)
   
   
   ##   selection of a balanced subset of edges for the assessment
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
      ##    computation of Balanced Error Rate
      BER.D2C<-BER(Ytrue,Yhat.D2C)
      BER.GS<-BER(Ytrue,Yhat.GS)
      BER.IAMB<-BER(Ytrue,Yhat.IAMB)
      
      cat("\n  DAG", r, ",", length(Ytrue), 
          " edges tested: BER.D2C=",BER.D2C, "BER.IAMB=",BER.IAMB,"BER.GS=",BER.GS,"\n")
   }
}



cat("\n Synthetic datasets: BER.D2C=",BER.D2C, "BER.IAMB=",BER.IAMB,"BER.GS=",BER.GS,"\n")




