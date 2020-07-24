
rm(list=ls())
library(D2C)
library(graph)
library(igraph)
library(gRbase)
library(ROCR)
library(RBGL)
library(bnlearn)


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
}

set.seed(0)

type="is.parent"
trainD2C<-NULL
knocked<-NULL

maxPar=4
wgt = 0.9

additive=FALSE
aBER=NULL
aBER2=NULL
aAUC=NULL
allDAG=NULL
for (iter in 0:100){
  nSamples=sample(50:200,1)
  if (iter %% 10 > 0 || iter==0){
    nNode=sample(5:20,1)
    g<-random_dag(1:nNode,maxpar=min(nNode,maxPar),wgt)
    cnt<-0
    while (sum(unlist(lapply(graph::edges(g),length)))<nNode & cnt<100){
      g<-random_dag(1:nNode,maxpar =min(nNode,maxPar),wgt)
      cnt<-cnt+1
      
    }
    G<-graph.adjacency(as(g,"matrix"))
    DAGg=as_graphnel(G)
  } else {
    DAGg=allDAG[which.max(aBER)]
  }
  allDAG=c(allDAG,DAGg)
  
  if (runif(1)<0.5){
    H = function() return(H_Rn(2)) #function() return(H_sigmoid(1))
  } else {
    H = function() return(H_Rn(1))
  }
  
  
  
  
  DAG = new("DAG.network",
            network=DAGg,H=H,additive=additive,
            weights=c(0.5,1),sdn=runif(1,0.2,0.5))
  
  observationsDAG <- compute(DAG,N=nSamples)
  cat("Observed data=",dim(observationsDAG), 
      "nEdges=",length(E(G)),"\n")
  
  trainDAG<-new("simulatedDAG",NDAG=1, N=100, noNodes=10,
                functionType = c("linear","quadratic","sigmoid"),
                seed=0,sdn=0.1,verbose=FALSE,
                additive=c(TRUE,FALSE),goParallel=FALSE)
  
  trainDAG@list.DAGs[[1]]=as_graphnel(G)
  trainDAG@list.observationsDAGs[[1]]=observationsDAG
  
  D2C<-new("D2C",sDAG=trainDAG,
           descr=descr,ratioEdges=0.2,
           max.features=20, type=type,goParallel=FALSE,
           verbose=FALSE,npar=1)
  
  
  
  
  if (iter>0){
    Nodes=nodes(DAGg)
    max.edges<-length(edgeList(DAGg))
    
    if (type=="is.parent"){
      subset.edges = matrix(unlist(edgeList(DAGg)),ncol=2,byrow = TRUE)
      subset.edges = unique(rbind(subset.edges,t(replicate(n =3*max.edges ,
                                                           sample(Nodes,size=2,replace = FALSE)))))
    } else {
      subset.edges = unique(t(replicate(n =4*max.edges ,sample(Nodes,size=2,replace = FALSE))))
    }
    
    
    Ahat.IAMB<-(amat(iamb(data.frame(observationsDAG),alpha=0.01,max.sx=3)))
    igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
    
    
    Yhat.D2C<-NULL
    Yhat.IAMB<-NULL
    phat.D2C<-NULL
    Ytrue<-NULL
    cat("D2C inferring", NROW(subset.edges), 
        "direct dependencies: please wait \n")
    
    for(jj in 1:NROW(subset.edges)){
      i=subset.edges[jj,1]
      j=subset.edges[jj,2]
      I =as(subset.edges[jj,1],"numeric")
      J =as(subset.edges[jj,2],"numeric") 
      if (length(intersect(c(I,J),knocked))==0){
        pred.D2C.rr =predict(trainD2C,I,J, observationsDAG,rep=4)$prob
        Yhat.D2C<-c(Yhat.D2C,round(pred.D2C.rr))
        phat.D2C<-c(phat.D2C,pred.D2C.rr)
        Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,I,J,type))
        Ytrue<-c(Ytrue,is.what(G,i,j,type)) 
        cat(".")
      }
    }
    
    
    aBER=c(aBER,BER(Ytrue,Yhat.D2C))
    aAUC=c(aAUC,AUC(Ytrue,phat.D2C))
    aBER2=c(aBER2,BER(Ytrue,Yhat.IAMB))
    cat("\n Assessment accuracy: \n BER D2C=",
        round(mean(aBER),2),"\n")
    cat("AUC D2C=",round(mean(aAUC),2),"\n")
    
    cat("\n  BER IAMB=",
        round(mean(aBER2),2),"\n")
    
    A=table(Ytrue,round(Yhat.D2C))
    rownames(A)=c("N","P")
    colnames(A)=c("N'","P'")
    print(A)
  }
  
  if (iter==0){
    allD2C<-D2C
  } else{
    allD2C@origX<-rbind(allD2C@origX,D2C@origX)
    allD2C@Y<-c(allD2C@Y,D2C@Y)
    cat("dim(allD2C@origX)=",dim(allD2C@origX),"NDAG computed=",iter,"\n")
  }
  
  
  trainD2C<-makeModel(allD2C,classifier="RF",EErep=2,verbose=FALSE)
  
  cat("D2C learner: \n # descriptors=",NCOL(allD2C@origX), 
      "\n # samples=",NROW(allD2C@origX), "\n # positives=",
      length(which(allD2C@Y==1)) , "\n")
  
  save(file="sequential.Rdata",list=c("allDAG","aBER","aAUC"))
}


