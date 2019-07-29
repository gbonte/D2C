rm(list=ls())
library("gbcode")
n=10
V=1:n
netwDAG<-new("graphNEL", nodes=as.character(1:n), edgemode="directed")
Nodes<-nodes(netwDAG)
NX<-Nodes[1]
NY<-Nodes[length(Nodes)]
netwDAG <- addEdge("2", NX, netwDAG, 1)
netwDAG <- addEdge(NX, NY, netwDAG, 1)
netwDAG <- addEdge("8", "2", netwDAG, 1)
netwDAG <- addEdge("3", NX, netwDAG, 1)
netwDAG <- addEdge("4", NY, netwDAG, 1)
netwDAG <- addEdge("5", "4", netwDAG, 1)
netwDAG <- addEdge("5", "3", netwDAG, 1)
netwDAG <- addEdge("9", "8", netwDAG, 1)
netwDAG <- addEdge("7", NY, netwDAG, 1)
netwDAG <- addEdge("7", NX, netwDAG, 1)
DAG = new("DAG.network",network=netwDAG,additive=TRUE,sdn=1)
N=100

nNodes <- numNodes(netwDAG)

topologicalOrder <-RBGL::tsort(netwDAG)

knockDAG<-netwDAG
wgt = runif(n = 1,min = 0.85,max = 1)
rDAG<-gRbase::random_dag(V,maxpar = 3,wgt) 
irDAG=igraph::graph.adjacency(as(rDAG,"matrix"))
plot(irDAG)

inetwDAG=igraph::graph.adjacency(as(netwDAG,"matrix"))
plot(inetwDAG)


setN<-NULL
for (nn in 2:(nNodes-1)){
  if ((! is.ancestor(inetwDAG,Nodes[nn],Nodes[1])) && ( is.ancestor(inetwDAG,Nodes[nn],Nodes[nNodes])) )
    setN<-c(setN,nn)
}
print(setN)
X=compute(DAG,N)
Yhatr=pred("rf",X[1:(N/2),c(1,setN)],X[1:(N/2),n],X[,c(1,setN)],class=FALSE)
Ytr=X[,n]
print(mean((Yhatr-Ytr)^2)/var(Ytr))

Xk=compute(DAG,100*N,knocked=1)
sX<-scale(rbind(X,Xk))

X<-sX[1:N,]
Xk<-sX[(N+1):NROW(sX),]



Yhat=pred("rf",X[,c(1,setN)],X[,n],Xk[,c(1,setN)],class=FALSE)
Y=Xk[,n]

Yhat2=pred("rf",X[,-n],X[,n],Xk[,-n],class=FALSE)


print(mean((Yhat-Y)^2)/var(Y))
print(mean((Yhat2-Y)^2)/var(Y))