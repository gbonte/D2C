
nNode=10
nSamples=100
g<-random_dag(1:nNode,maxpar=3,wgt=1)


G<<-graph.adjacency(as(g,"matrix"))


if (runif(1)<0.5){
  H = function() return(H_sigmoid(1))
} else {
  H = function() return(H_Rn(1))
}

additive=sample(c(TRUE,FALSE),1)

DAG = new("DAG.network",
          network=as_graphnel(G),H=H,additive=additive,sdn=runif(1,0.2,0.5))

DX<-NULL
DY<-NULL
knockset<-c(1,3,7)
for (k in knockset){

  X <- compute(DAG,N=nSamples,knocked=k)
  
  for (i in setdiff(nNode,k))
    for (j in setdiff(nNode,k)){
      DX<-rbind(DX,c(npred(X[,k],X[,i]),npred(X[,k],X[,j]),norminf(X[,k],X[,i],X[,j]), norminf(X[,k],X[,j],X[,i])))
      if (is.parent(i,j))
        DY<-c(DY,1)
      else
        DY<-c(DY,0)
      
      
    }


}

F <- randomForest(x =X ,y = factor(Y),importance=TRUE)

F
