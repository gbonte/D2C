rm(list=ls())

library(graph)
library(igraph)
library(readxl)
library(gRbase)
library(D2C)
library(ROCR)
library(RBGL)


DX<-NULL
DY<-NULL
lin<-TRUE
for (r in 1:30){
  set.seed(r)
  nNode=sample(21:50,1)
  nSamples=10
  g<-random_dag(1:nNode,maxpar=3,wgt=0.9)
  G<-graph.adjacency(as(g,"matrix"))
  
  if (runif(1)<0.5){
    H = function() return(H_sigmoid(1))
  } else {
    H = function() return(H_Rn(1))
  }
  
  additive=sample(c(TRUE,FALSE),1)
  
  DAG = new("DAG.network",
            network=as_graphnel(G),H=H,additive=additive,sdn=runif(1,0.2,0.5))
  
  
  knockset<-sample(1:nNode,10)
  X <- compute(DAG,N=nSamples,knock=knockset)
  for (k in knockset){
    
    
    
    for (i in sample(setdiff(1:nNode,knockset),5))
      for (j in sample(setdiff(1:nNode,knockset),5)){
        DX<-rbind(DX,c(descriptor(X,i,j,ns=2,biv=FALSE),c(npred(X[,k],X[,i]),npred(X[,k],X[,j]),norminf(X[,k],X[,i],X[,j]), norminf(X[,k],X[,j],X[,i]))))
       ##DX<-rbind(DX,c(npred(X[,i],X[,j],lin=lin),npred(X[,k],X[,i],lin=lin),npred(X[,k],X[,j],lin=lin),
         ##              norminf(X[,k],X[,i],X[,j],lin=lin), norminf(X[,k],X[,j],X[,i],lin=lin)))
        
        if (is.parent(G,as.character(i),as.character(j)))
          DY<-c(DY,1)
        else
          DY<-c(DY,0)
        
        cat(".")
      }
    
    
  }
  print(r)
}

w0<-which(DY==0)
w1<-which(DY==1)
if (length(w0)>length(w1))
  w0<-sample(w0,length(w1))

if (length(w1)>length(w0))
  w1<-sample(w1,length(w0))
wc=which(apply(DX,2,sd)<0.01)
if (length(wc)>0)
  DX=DX[,-wc]
DX<-scale(DX[c(w0,w1),])
DY<-DY[c(w0,w1)]

print(dim(DX))
F <- randomForest(x =DX ,y = factor(DY),importance=TRUE)

F
