library(D2C)
## script showing the link between multivariate series (3 series) and DAG


set.seed(0)

G<-genSTAR(n=3,nn=7,NN=100,sd=0.5,num=1,loc=1,verbose=TRUE)
## n=3 series represented with nn=5 nodes each
## data are sampled by the generator num
## loc is the size of the neighborhood 

## NODES ENCODING (3*7 nodes)
## nodes 1:nn= series 1  1:y1[t]   2:y1[t-1] 3:y1[t-2] 4:y1[t-3] 5:y1[t-4] 6:y1[t-5 ] 7:y1[t-6 ]
## nodes (nn+1):(2*nn)= series 2: 8:y2[t]  9:y2[t-1]  10:y2[t-2]  11:y2[t-3] 12:y2[t-4] 13:y2[t-5] 14:y2[t-6]
## nodes (2*nn+1):(3*nn)= series 3: 15:y3[t]  16:y3[t-1]  17:y3[t-2]  18:y3[t-3]  19:y3[t-4] 20:y3[t-5] 21:y3[t-6]

print(G$fs) ## number of parents for each node in G$doNeigh
# 0 2

print(G$doNeigh)
                       
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
