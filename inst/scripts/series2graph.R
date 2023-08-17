library(D2C)
library(igraph)
## script showing the link between multivariate series and DAG

set.seed(0)
G<-genSTAR(n=2,nn=5,NN=100,sd=0.5,num=1,loc=1,verbose=TRUE)
## n=2 series represented with nn=5 nodes each
## data are sampled by the generator num=1
## loc is the size of the neighborhood 

## nodes 1:nn= series 1  1:y1[t] 2:y1[t-1] 3:y1[t-2] 4:y1[t-3] 5:y1[t-4]
## nodes (nn+1):(2*nn)= series 2: 6:y2[t]  7:y2[t-1]  8:y2[t-2]  9:y2[t-3]    10:y2[t-4]

print(G$fs) ## lag of parents for each node in G$doNeigh
# 0 3
print(G$doNeigh) # same size as n
## G$doNeigh[[1]] denotes to which series belong the parents of the 1st series
#[[1]]
#[1] 1 : means +1 so the parents are in the series y(1+1)=y2
#           y1[t] (N1) depends on y2[t-fs[1]-1]=y2[t-1] (N7)  and y2[t-fs[2]-1]=y2[t-4] (N10)
#               N1 <- N7               N1<-N10

## G$doNeigh[[2]] denotes to which series belong the antecedent of the 2nd series
#[[2]]
#[1] 0  : the parents are in the same series y2
#       y2[t] (N6) depends on y2[t-fs[1]-1]=y2[t-1] (N7)  and y2[t-fs[2]-1]=y2[t-4] (N10)
# -1  : y2[t] (N6) depends on y1[t-fs[1]-1]=y1[t-1] (N2)  and y1[t-fs[2]-1]=y1[t-4] (5)
##        N6 <- N7  N6 <-10 N6 <-N2  N6<-N5

## 15 edges
## y1[t]=f(y2[t-1],y2[t-4])   N1<-N7 N1 <-N10
## y1[t-1]=f(y2[t-2],y2[t-5]) N2<-N8
## y1[t-2]=f(y2[t-3],y2[t-6]) N3<-N9
## y1[t-3]=f(y2[t-4],y2[t-7]) N4<-N10
## y1[t-4]=f(y2[t-5],y2[t-8]) 
## y2[t]=f(y1[t-1],y1[t-4],y2[t-1],y2[t-4]) N6<-N2, N6<-N5 N6<-N7 N6<-N10
## y2[t-1]=f(y1[t-2],y1[t-5],y2[t-2],y2[t-5]) N7<-N3 N7<-N8 
## y2[t-2]=f(y1[t-3],y1[t-6],y2[t-3],y2[t-6]) N8<-N4 N8<-N9
## y2[t-3]=f(y1[t-4],y1[t-7],y2[t-4],y2[t-7]) N9<-N5 N9<-N10
## y2[t-4]=f(y1[t-5],y1[t-7],y2[t-4],y2[t-7])

print(G$DAG)

## number of lags 
print(dim(G$D))

igraphDAG<-graph.adjacency(as(G$DAG,"matrix"))
plot(igraphDAG,layout=layout_with_sugiyama)