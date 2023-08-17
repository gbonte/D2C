library(D2C)
## script showing the link between multivariate series (3 series) and DAG


set.seed(0)

G<-genSTAR(n=3,nn=5,NN=100,sd=0.5,num=1,loc=1,verbose=TRUE)
## n=3 series represented with nn=5 nodes each
## data are sampled by the generator num
## loc is the size of the neighborhood 

## NODES ENCODING (3*15 nodes)
## nodes 1:nn= series 1  1:y1[t]   2:y1[t-1] 3:y1[t-2] 4:y1[t-3] 5:y1[t-4]
## nodes (nn+1):(2*nn)= series 2: 6:y2[t]  7:y2[t-1]  8:y2[t-2]  9:y2[t-3]    10:y2[t-4]
## nodes (2*nn+1):(3*nn)= series 3: 11:y3[t]  12:y3[t-1]  13:y3[t-2]  14:y3[t-3]    15:y3[t-4]

print(G$fs) ## number of parents for each node in G$doNeigh
# 0 2

print(G$doNeigh)
## G$doNeigh[[1]] denotes to which series belong the parents of the 1st series
#[[1]]
#[1] 0 y1[t] (N1) depends on y1[t-fs[1]-1]=y1[t-1] (N2)  and y1[t-fs[2]-1]=y1[t-3] (N4)
#               N1 <- N2               N1<-N4

## G$doNeigh[[2]] denotes to which series belong the antecedent of the 2nd series
#[[2]]
#[1] 1  : y2[t] (N6) depends on y3[t-fs[1]-1]=y3[t-1] (N12)  and y3[t-fs[2]-1]=y3[t-3] (N14)
# -1  : y2[t] (N6) depends on y1[t-fs[1]-1]=y1[t-1] (N2)  and y1[t-fs[2]-1]=y1[t-3] (N4)
##        N6 <- N12  N6 <-14 N6 <-N2  N6<-N4

## G$doNeigh[[3]] denotes to which series belong the antecedent of the 3rd series
#[[2]]
#[1] 0  : y3[t] (N11) depends on y3[t-fs[1]-1]=y3[t-1] (N12)  and y3[t-fs[2]-1]=y3[t-3] (N14)
##        N11 <- N12  N11 <-14 


## 24 Edges
## y1[t]=f(y1[t-1],y1[t-3]) N1<-N2, N1<-N4
##                          N2<-N3 N2<-N5
##                          N3<-N4
##                          N4<-N5
## y2[t]=f(y3[t-1],y3[t-3],y1[t-1],y1[t-3])  N6<-N12  N6<-N14  N6<-N2  N6<-N4
##                                        N7<-N13  N7<-N15  N7<-N3  N7<-N5
##                                        N8<-N14  N8<-N4 
##                                        N9<-N15  N9<-N5
## y3[t]=f(y3[t-1],y3[t-3])     N11<-N12  N11<-N14
##                              N12<-N13  N12<-N15
##                              N13<-N14  
#                               N14<-N15
print(G$DAG)

## number of lags 
print(dim(G$D))

igraphDAG<-graph.adjacency(as(G$DAG,"matrix"))
plot(igraphDAG,layout=layout_with_sugiyama)


#layout_as_star,layout_as_tree, layout_in_circle(), layout_nicely(),
#  layout_on_grid(), layout_on_sphere(), layout_randomly(), 
# layout_with_dh(), layout_with_fr(), layout_with_gem(), 
## layout_with_graphopt(), 
## layout_with_kk(), layout_with_lgl(), layout_with_mds(), 
## layout_with_sugiyama()
