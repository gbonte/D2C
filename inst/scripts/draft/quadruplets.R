#library(D2C)
#library(gRbase)
#library(Rgraphviz)
#library(igraph)

graphics.off()


noNodes=3
V=1:noNodes
maxpar = 3
wgt = runif(n = 1,min = 0.85,max = 1)
H = function() return(H_Rn(1))
netwDAG<-random_dag(V,maxpar = maxpar,wgt)  
nodes(netwDAG)<-as.character(V)

### random_dag {gRbase}: generate a graphNEL random directed acyclic graph (DAG)

netwDAG<-random_dag(V,maxpar = 3,1)


plot(netwDAG)


iDAG=graph.adjacency(as(netwDAG,"matrix")) 

DAG = new("DAG.network",
          network=netwDAG,H=H,additive=TRUE,sdn=0.1)


X = compute(DAG,N=100)

n<-NCOL(X)

Nodes=V(iDAG)
for (i in 1:(n-1))
  for (j in (i+1):n){
    
    cat(Nodes[i], "is.parent", Nodes[j], "?:",is.parent(iDAG,Nodes[i],Nodes[j]),"\n")
    cat(Nodes[j], "is.parent", Nodes[i], "?:",is.parent(iDAG,Nodes[j],Nodes[i]),"\n")
  }
