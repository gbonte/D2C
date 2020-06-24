rm(list=ls())
library(D2C)
type="is.parent"

is.what<-function(iDAG,i,j){
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
library(vars)
namefile<-"../D2Cdata2/traintestSTAR.RData"
load(namefile)
BERD2C<-NULL
BERN<-NULL
BERN2<-NULL
Strue<-NULL
Shat<-NULL
Shat2<-NULL

for (r in 1:200){
  observedData<-testDAG@list.observationsDAGs[[r]]
  data<-testDAG@list.Y[[r]]
  fs<-testDAG@list.fs[[r]]
  Y <- data.frame(scale(data))
  
  n <- ncol(data)
  if (n<=5){
    print(n)
    doNeigh=testDAG@list.doNeigh[[r]]
    print(doNeigh)
    trueDAG<-testDAG@list.DAGs[[r]]
    graphTRUE<- gRbase::as.adjMAT(trueDAG)
    igraph.TRUE<-igraph::graph.adjacency(graphTRUE[as.character(1:NCOL(graphTRUE)),as.character(1:NCOL(graphTRUE))])
    
    plot(igraph.TRUE)
    
    nn=6
    M=MakeEmbedded(Y,n=numeric(n)+nn,delay=numeric(n),hor=rep(1,n),w=1:n)
    
    D=M$inp
    P<-array(0,c(n,n*nn))
    PP<-array(0,c(n,n*nn))
    ## P[i,j] probability that i is caused by j
    for (i in seq(1,NCOL(D),by=nn)){
      for (j in (1:NCOL(D))){
        if (j%%nn != 1){
          cat(".")
          #i=",i,"j=",j,"\n")
          
          P[ceiling(i/nn),j]=predict(trainD2C,j,i, D,rep=3)$prob
          PP[ceiling(i/nn),j]=is.what(igraph.TRUE,j,i)
        }
        
      }
    }
    
    BERD2C<-c(BERD2C,BER(c(PP),c(round(P))))
    cat("\n Err=",mean(BERD2C),"\n")
    
    ## CauseMe requires to upload a score matrix and
    ## optionally a matrix of p-values and time lags where
    ## the links occur
    
    ## In val_matrix an entry [i, j] denotes the score for the link i --> j and
    ## must be a non-negative real number with higher values denoting a higher
    ## confidence for a link.
    ## Fitting a VAR model results in several lagged coefficients for a
    ## dependency of j on i.
    ## Here we pick the absolute value of the coefficient corresponding to the
    ## lag with the smallest p-value.
    
    val_matrix <- array(rep(0, n*n), c(n, n))
    V <- array(rep(0, n*n), c(n, n))
    ## Matrix of p-values
    p_matrix <-  array(rep(0, n*n), c(n, n))
    
    ## Matrix of time lags
    lag_matrix <-  array(rep(0, n*n), c(n, n))
    loc=1
    for (i in 1:n){
      neigh=i+doNeigh[[i]]
      V[neigh,i]=1
      for (j in 1:n){
        
        val_matrix[i, j] <- max(P[j,((i-1)*nn+1):(i*nn)])
        ## p_matrix[i, j] <- pvalues[i, j, min_pvalue_lag]  # Store the p-value
        lag_matrix[i, j] <- which.max(P[j,((i-1)*nn+1):(i*nn)])  # Store the time lag 
      }
    }
    
    
    
    scores <-  as.vector(val_matrix)
    print(val_matrix)
    Strue<-c(Strue,c(V))
    Shat<-c(Shat,c((val_matrix)))
    BERN<-c(BERN,BER(c(V),c(round(val_matrix))))
    cat("BErr N=",mean(BERN),"AUC=",AUC(Strue,Shat), "\n ------- \n" )
    ################## VAR method
    maxlags=5
    N <- ncol(data)
    results <- VAR(data, p = maxlags)  # Fit the VAR model
    values = array(rep(NaN, N*N*maxlags), c(N, N, maxlags))  # Store the coeff
    pvalues = array(rep(NaN, N*N*maxlags), c(N, N, maxlags))  # Store the p-values
    
    for (i in 1:N){  # Fill the value and pvalue matrix
      values[i, , ] <- array(results$varresult[[i]]$coefficients[1:(N*maxlags)], c(N, maxlags))
      
      pvalues[i, ,] <- array(coef(results)[[i]][, 4], c( c(N, maxlags)))
    }
    
    # CauseMe requires to upload a score matrix and
    # optionally a matrix of p-values and time lags where
    # the links occur
    
    # In val_matrix an entry [i, j] denotes the score for the link i --> j and
    # must be a non-negative real number with higher values denoting a higher
    # confidence for a link.
    # Fitting a VAR model results in several lagged coefficients for a
    # dependency of j on i.
    # Here we pick the absolute value of the coefficient corresponding to the
    # lag with the smallest p-value.
    
    val_matrix <- array(rep(NaN, N*N), c(N, N))
    
    # Matrix of p-values
    p_matrix <-  array(rep(NaN, N*N), c(N, N))
    
    # Matrix of time lags
    lag_matrix <-  array(rep(NaN, N*N), c(N, N))
    
    for (i in 1:N){
      for (j in 1:N){
        min_pvalue_lag <- match(min(pvalues[i,j,]), pvalues[i,j,])  # Find the min p_value
        val_matrix[i, j] <- abs(values[i, j, min_pvalue_lag])  # Set the actual absolute value at that time-lag
        p_matrix[i, j] <- pvalues[i, j, min_pvalue_lag]  # Store the p-value
        lag_matrix[i, j] <- min_pvalue_lag  # Store the time lag 
      }
    }
    Shat2<-c(Shat2,c((val_matrix)))
    BERN2<-c(BERN2,BER(c(V),c(round(val_matrix))))
    cat("BErr VAR N=",mean(BERN2),"AUC=",AUC(Strue,Shat2),"\n ------- \n" )
  }
}
## The A_ij element of this matrix, corresponding to the i-th row and j-th column, 
##indicates the score of a causal link i to j, which is a non-negative real number that can either be 
## a probability estimates, confidence value, or binary decision. Higher values indicate more confidence in a link. 
## If the absence of a link is estimated with maximum certainty, then A_ij = 0. Real values between 0 and 1 indicate 
## probabilities in between. The matrix must be stored flattened in row-major (C-style) order. T

## pvalues and lags are recommended for a more comprehensive method evaluation,
## but not required. Then you can leave the dictionary field empty  



