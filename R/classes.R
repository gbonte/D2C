#' @import RBGL gRbase randomForest xgboost Rgraphviz methods foreach kernlab MASS igraph graph e1071
 

 

#########################################
########   class D2C.descriptor
#########################################

##' An S4 class to store the descriptor parameters
setClass("D2C.descriptor",
         slots = list(lin="logical", acc="logical",
                      struct="logical",pq="numeric",
                      bivariate="logical",residual="logical",
                      stabD="logical",
                      diff="logical",ns="numeric",
                      maxs="numeric",boot="character"))

##' creation of a D2C.descriptor
##' @name D2C descriptor
##' @param .Object : the D2C.descriptor object
##' @param lin \{TRUE, FALSE\}: if TRUE it uses a linear model to assess a dependency, otherwise a local learning algorithm (lazy package)
##' @param acc \{TRUE, FALSE\}: if TRUE it uses the accuracy of the regression as a descriptor
##' @param struct	\{TRUE, FALSE\}: if TRUE it uses the ranking in the Markov Blanket as a descriptor
##' @param pq :a vector of quantiles used to compute the descriptors
##' @param bivariate \{TRUE, FALSE\}: if TRUE it includes also the descriptors of the bivariate dependence
##' @param residual \{TRUE, FALSE\}: if TRUE it includes also the residual in the descriptor computation
##' @param diff \{TRUE, FALSE\}: if TRUE it includes also the difference values in the descriptor computation (only for time series)
##' @param stabD \{TRUE, FALSE\}: if TRUE it includes also the stability in the descriptor computation 
##' @param ns : size of the Markov Blanket returned by the mIMR algorithm
##' @param boot : bootstrap algorithm
##' @param maxs : max size of distribution samples
##' @references Gianluca Bontempi, Maxime Flauder (2015) From dependency to causality: a machine learning approach. JMLR, 2015, \url{http://jmlr.org/papers/v16/bontempi15a.html}
##' @examples
##' require(RBGL)
##' require(gRbase)
##' require(foreach)
##'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
##'
##' @export
setMethod("initialize",
          "D2C.descriptor",
          function(.Object, lin=TRUE, acc=TRUE,
                   struct=FALSE,pq=c(0.1, 0.25, 0.5, 0.75, 0.9),
                   bivariate=FALSE,ns=4,boot="rank",maxs=20,diff=FALSE, residual=FALSE, 
                   stabD=FALSE)
          {
            
            .Object@lin <- lin
            .Object@acc <- acc
            .Object@struct <- struct
            .Object@bivariate <- bivariate
            .Object@pq <- pq
            .Object@ns <- ns
            .Object@maxs <- maxs
            .Object@boot <- boot
            .Object@residual <- residual
            .Object@diff <- diff
            .Object@stabD <- stabD
            .Object
          }
)




#########################################
########   class DAG.network
#########################################




##' An S4 class to store DAG.network
##' @param network : object of class "graph"
setClass("DAG.network",  slots = list(network = "graph",additive="logical",maxV="numeric",
                                      exosdn="numeric"))


##' creation of a DAG.network
##' @name DAG network
##' @param network : object of class "igraph"
##' @param sdn : error noise
##' @param exosdn :  sdvar  exogenous inputs
##' @param sigma : function returning the additive noise
##' @param H : function describing the type of the dependency.
##' @param additive : if TRUE the output is the sum of the H transformation of the inputs, otherwise it is the H transformation of the sum of the inputs.
##' @param weights : [lower,upper], lower and upper bound of the values of the linear weights
##' @param maxV: maxim absolute value
##' @export
setMethod("initialize", signature="DAG.network",
          function(.Object, network,
                   sdn=0.5,
                   exosdn=1,
                   sigma=function(x) {
                     return(rnorm(n = 1,sd = sdn))},
                   H=c(function(x) return(H_Rn(1)),
                       function(x) return(H_Rn(2))),
                   additive= TRUE,
                   weights=c(0.8,2),maxV=5){
            DAG = network
            .Object@additive=additive
            .Object@maxV=maxV
            .Object@exosdn=exosdn
            if(!is.DAG(DAG)) {
              stop("it is not a DAG")
            } else  {
              nodeDataDefaults(DAG,"bias") <-0
              nodeDataDefaults(DAG,"sigma") <-sigma
              nodeDataDefaults(DAG,"seed") <-NA
              edgeDataDefaults(DAG,"H") <- function(x) return(x)
              for( n in nodes(DAG)){
                set.seed(as.numeric(Sys.time()))
                nodeData(DAG,n=n,"seed")<-runif(1,1,10000)
                nodeData(DAG,n=n,"sigma")<-function(x) {
                  return(rnorm(n = 1,sd = runif(1,0.9*sdn,sdn)))}
              }
              for( edge in edgeList(DAG)){
                edgeData(DAG, from=edge[1], to=edge[2], attr="weight") <- runif(1,weights[1],weights[2])*sample(c(-1,1),1)
                ## setting of random linear weights within the specified bounds
                
                if (length(H)>1){
                  Hi<-sample(H,1)[[1]]
                }else{
                  Hi<-H
                }
                
                edgeData(DAG, from=edge[1], to=edge[2],attr="H") <- Hi()
                
              }
            }
            .Object@network <- DAG
            
            return(.Object)
          }
)




#' @docType methods
setGeneric("compute", function(object,...) {standardGeneric("compute")})
##' generate N samples according to the network distribution
##' @name compute
##' @param N: the number of samples generated according to the network
##' @param object: a DAG.network object
##' @return a N*nNodes matrix
##' @export
setMethod("compute", signature="DAG.network",  
          function(object, N=50){
            if(!is.numeric(N))
              stop("N is not numeric")
            save.seed <- get(".Random.seed", .GlobalEnv)
            DAG = object@network
            maxV= object@maxV
            nNodes <- numNodes(DAG)
            
            topologicalOrder <-tsort(DAG)
            
            
            DD<-NULL
            Nsamples<-0
            it<-0
            while (Nsamples < N & it <2){
              D <- matrix(NA,nrow=2*N,ncol=nNodes)
              colnames(D) <- 1:nNodes
              it<-it+1
              for (i in topologicalOrder){
                bias = nodeData(DAG,n=i,attr="bias")[[1]]
                sigma = nodeData(DAG,n=i,attr="sigma")[[1]]
                seed = nodeData(DAG,n=i,attr="seed")[[1]]+it
                inEdg <-  inEdges(node=i,object=DAG)[[1]]
                
                if (length(inEdg)==0 ){
                  set.seed(seed)
                  dr=rnorm(2*N,sd=object@exosdn)
                  D[,i]<-bias + dr  #replicate(N,sigma())
                } else  {
                  D[,i]<-bias
                  Xin<-NULL
                  for (j in  inEdg){
                    ## it computes the linear combination of the inputs
                    inputWeight = edgeData(self=DAG,from=j,to=i,attr="weight")[[1]]
                    H = edgeData(self=DAG,from=j,to=i,attr="H")[[1]]
                    
                    if (object@additive){
                      D[,i]<- D[,i] + H(D[,j]) *  inputWeight
                    }else{
                      ## it stacks inputs in a matrix
                      Xin<-cbind(Xin,D[,j]*  inputWeight)
                      
                    }
                  }
                  if (!object@additive){
                    H = edgeData(self=DAG,from=inEdg[1],to=i,attr="H")[[1]]
                    if (length(inEdg)==1)
                      D[,i]<-  H(Xin)
                    else
                      D[,i]<-  H(apply(Xin,1,sum))
                    
                  }
                  set.seed(seed)
                  
                  D[,i] <- (D[,i] + replicate(2*N,sigma())/it)  ## additive random noise
                  
                }
              } ## for i
              #col.numeric<-as(colnames(D),"numeric")
              #D<-D[,topologicalOrder[order(col.numeric)]]
              Dmax<-apply(abs(D),1,max)
              wtoo<-union(which(Dmax>maxV),which(is.na(Dmax)))
              if (length(wtoo)>0)
                D=D[-wtoo,]
              
              DD<-rbind(DD,D) ## remove divergent samples
              Nsamples=NROW(DD)
              #print(Nsamples) 
              #if (Nsamples==0)
              #  browser()
            }
            assign(".Random.seed", save.seed, .GlobalEnv)
            
            if (N==0)
              return(DD)
            if (N<=NROW(DD))
              DD=DD[1:N,]
            return(DD)
          })

#' @docType methods
setGeneric("counterfact", function(object,...) {standardGeneric("counterfact")})
##' generate N samples according to the network distribution by modifying the original dataset
##' @name counterfact
##' @param DN: original dataset
##' @param knocked: the set of manipulated (e.g. knocked genes) nodes 
##' @param object: a DAG.network object
##' @return a N*nNodes matrix
##' @export
setMethod("counterfact", signature="DAG.network",  
          function(object, DN, delta,
                   knocked=NULL){
            if(!is.numeric(N))
              stop("N is not numeric")
            save.seed <- get(".Random.seed", .GlobalEnv)
            DAG = object@network
            nNodes <- numNodes(DAG)
            maxV= object@maxV
            topologicalOrder <-tsort(DAG)
            posknock<-min(match(knocked,topologicalOrder))
            beforeknock<-numeric(nNodes)
            if (posknock>1)
              beforeknock[1:(posknock-1)]=1
            D <- DN
            N<-NROW(DN)
            
            
            for (ii in 1:length(topologicalOrder)){
              
              i=topologicalOrder[ii]
              if (beforeknock[ii]==0){
                
                bias = nodeData(DAG,n=i,attr="bias")[[1]]
                sigma = nodeData(DAG,n=i,attr="sigma")[[1]]
                seed = nodeData(DAG,n=i,attr="seed")[[1]]
                inEdg <-  inEdges(node=i,object=DAG)[[1]]
                
                if (length(inEdg)==0 || is.element(i,knocked)){
                  if (length(inEdg)==0){
                    D[,i]<-DN[,i] 
                  }
                  if (is.element(i,knocked))
                    D[,i]<-delta
                } else  {  ##  if (length(inEdg)==0
                  D[,i]<-bias
                  Xin<-NULL
                  for(j in  inEdg)
                    ## it computes the linear combination of the inputs
                  {
                    inputWeight = edgeData(self=DAG,from=j,to=i,attr="weight")[[1]]
                    H = edgeData(self=DAG,from=j,to=i,attr="H")[[1]]
                    
                    if (object@additive){
                      D[,i]<- D[,i] + H(D[,j]) *  inputWeight
                    }else{
                      Xin<-cbind(Xin,D[,j]*  inputWeight)
                      ##D[,i]<- D[,i] + (D[,j]) *  inputWeight
                    }
                  }
                  if (!object@additive){
                    
                    H = edgeData(self=DAG,from=inEdg[1],to=i,attr="H")[[1]]
                    if (length(inEdg)==1)
                      D[,i]<-  H(Xin)
                    else
                      D[,i]<-  H(apply(Xin,1,sum))
                    
                  }
                  set.seed(seed) 
                  
                  D[,i] <- D[,i] + replicate(N,sigma())  ## use of sigmoid function to saturate + additive random noise
                  
                }
              } # if beforeknock
            } 
            
            assign(".Random.seed", save.seed, .GlobalEnv)
            return(D)
            
          })


#########################################
########   class simulatedDAG
#########################################

##' An S4 class to store a list of DAGs and associated observations
##' @name simulatedDAG
##' @param list.DAGs : list of stored DAGs
##' @param list.observationsDAGs : list of observed datasets, each sampled from the corresponding member of list.DAGs
##' @param NDAG  : number of DAGs.
##' @param functionType : type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid")
##' @param seed : random seed
setClass("simulatedDAG",
         slots = list(list.DAGs="list",list.observationsDAGs="list",
                      NDAG="numeric", functionType="character",seed="numeric",
                      additive="logical",weights="numeric"))




##' creation of a "simulatedDAG" containing a list of DAGs and associated observations
##' @param .Object : simulatedDAG object
##' @param NDAG : number of DAGs to be created and simulated
#' @param noNodes  : number of Nodes of the DAGs. If it is a two-valued vector , the value of Nodes is randomly sampled in the interval
#' @param N  : number of sampled observations for each DAG. If it is a two-valued vector [a,b], the value of N is randomly sampled in the interval [a,b]
#' @param sdn : standard deviation of aditive noise. If it is a two-valued vector, the value of N is randomly sampled in the interval
#' @param seed : random seed
#' @param verbose : if TRUE it prints out the state of progress
#' @param functionType : type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid","kernel")
#'  @param quantize  : if TRUE it discretize the observations into two bins. If it is a two-valued vector [a,b], the value of quantize is randomly sampled in the interval [a,b]
#'  @param maxpar.pc  : maximum number of parents expressed as a percentage of the number of nodes
#'  @param goParallel : if TRUE it uses parallelism
#'  @param additive : if TRUE the output is the sum of the H transformation of the inputs, othervise it is the H transformation of the sum of the inputs.
#' @param weights : [lower,upper], lower and upper bound of the values of the linear weights
#'  @param maxV : max accepted value
#'  @references Gianluca Bontempi, Maxime Flauder (2015) From dependency to causality: a machine learning approach. JMLR, 2015, \url{http://jmlr.org/papers/v16/bontempi15a.html}
#' @examples
#' require(RBGL)
#' require(gRbase)
#' require(foreach)
#' descr=new("D2C.descriptor")
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedDAG",NDAG=10, N=c(50,100),noNodes=c(15,40),
#'              functionType = "linear", seed=0,sdn=c(0.45,0.75))
#' @export
#'
#'
setMethod("initialize",
          "simulatedDAG",
          function(.Object, NDAG=1,
                   noNodes=sample(10:20,size=1),functionType="linear",
                   quantize=FALSE,maxpar.pc=0.05,
                   verbose=TRUE,N=sample(100:500,size=1),
                   seed=1234,sdn=0.5, goParallel=FALSE,additive=FALSE,
                   weights=c(0.5,1),maxV=5)
          {
            
            ##generate a training set
            ## NDAG the number of network to use
            ##functionType example : "R1" "R2" "sigmoid1"
            
            
            if (goParallel){
              # cl <- makeForkCluster(10)
              #  registerDoParallel(cl)
              `%op%` <-  `%dopar%` 
            }   else {
              `%op%` <-`%do%`
            }
            .Object@functionType=functionType
            .Object@seed=seed
            .Object@weights=weights
            X=NULL
            Y=NULL
            list.DAGs=NULL
            list.observationsDAGs=NULL
            if (NDAG<=0)
              return(.Object)
            
            
            
            FF<-foreach (i=1:NDAG) %op%{
              ##for (i in 1:NDAG){
              set.seed(seed+i)
              
              N.i<-N
              if (length(N)>1)
                N.i<-sample(N[1]:N[2],1)
              
              quantize.i<-quantize
              if (length(quantize)>1)
                quantize.i<-sample(quantize,1)
              
              noNodes.i<-max(3,noNodes[1])
              if (length(noNodes)==2)
                noNodes.i<-sample(max(3:noNodes[1]):max(3:noNodes[2]),1)
              
              
              sdn.i<-sdn
              if (length(sdn)>1)
                sdn.i<-runif(1,sdn[1],sdn[2])
              
              functionType.i<-functionType
              if (length(functionType.i)>1)
                functionType.i<-sample(functionType,1)
              
              additive.i<-additive
              if (length(additive)>1)
                additive.i<-sample(additive,1)
              
              
              sdn.ii=sdn.i
              weights.i=weights
              maxV.i=maxV
              
              HH<-NULL
              for (functionType.i in functionType){
                if(functionType.i=="linear"){
                  H = function() return(H_Rn(1))
                  
                }else if(functionType.i=="quadratic"){
                  H = function() return(H_Rn(2))
                  
                }else if(functionType.i=="sigmoid"){
                  H = function() return(H_sigmoid(1))
                  
                } else if(functionType.i=="kernel"){
                  H = function() return(H_kernel())
                  
                }
                HH<-c(HH,H)
              }
              
              cnt2=0
              while(1){
                V=1:max(4,noNodes.i-cnt2)
                
                maxpar.pc.i<-pmin(0.99,maxpar.pc)
                if (length(maxpar.pc.i)>1)
                  maxpar.pc.i<-sample(maxpar.pc,1)
                
                maxpar = round(maxpar.pc.i*noNodes)
                
                
                wgt = runif(n = 1,min = 0.85,max = 1)
                
                
                netwDAG<-random_dag(V,maxpar = maxpar,wgt)  
                ### random_dag {gRbase}: generate a graphNEL random directed acyclic graph (DAG)
                
                nodes(netwDAG)<-as.character(V)
                
                cnt<-1
                while (sum(unlist(lapply(graph::edges(netwDAG),length)))<2 ){
                  maxpar = sample(1:max(3,round(noNodes.i/3)),size=1)
                  netwDAG<-random_dag(V,maxpar = maxpar,1)
                  nodes(netwDAG)<-as.character(V)
                  cnt<-cnt+1
                  if (cnt>50){
                    netwDAG<-new("graphNEL", nodes=as.character(1:noNodes.i), edgemode="directed")
                    netwDAG <- addEdge("1", "3", netwDAG, 1)
                    netwDAG <- addEdge("2", "3", netwDAG, 1)
                  }
                }
                
                ## it iterates until there is a dataset sufficiently large
                
                
                
                set.seed(as.numeric(Sys.time()))
                
                DAG = new("DAG.network",
                          network=netwDAG,H=HH,additive=additive.i,
                          sdn=sdn.ii,weights=weights.i,maxV=max(maxV.i,1))
                HH = function() return(H_Rn(1))
                sdn.ii=0.99*sdn.ii
                weights.i=0.99*weights.i
                maxV.i=maxV.i-1
                cnt2=cnt2+1
                
                observationsDAG = compute(DAG,N=max(20,N.i-cnt2))
                if (! (is.null(observationsDAG) | is.vector(observationsDAG)))
                  if (NROW(observationsDAG)>max(10,round(N.i/2)-cnt2))
                    break;
              }
              
              if (quantize.i)
                observationsDAG<-apply(observationsDAG,2,quantization)
              
              if (any(is.na(observationsDAG) ))
                stop("simulatedDAG: NA in data generation")
              if (any(is.infinite(observationsDAG)))
                stop("simulatedDAG: inf in data generation")
              if (verbose){
                
                cat("simulatedDAG: DAG number:",i,"generated: #nodes=", length(V),
                    "# edges=",sum(unlist(lapply(graph::edges(netwDAG),length))), "# samples=", NROW(observationsDAG), "\n")
                
              }
              
              
              
              list(observationsDAG=observationsDAG,netwDAG=netwDAG)
              
            } ## foreach
            
            
            .Object@list.DAGs=lapply(FF,"[[",2)
            .Object@list.observationsDAGs=lapply(FF,"[[",1)
            
            to.remove=which(unlist(lapply(lapply(.Object@list.DAGs,edgeList),length))==0)
            if (length(to.remove)>0){
              .Object@list.DAGs=.Object@list.DAGs[-to.remove]
              .Object@list.observationsDAGs=.Object@list.observationsDAGs[-to.remove]
            }
            .Object@NDAG=length(.Object@list.DAGs)
            
            .Object
          }
)



#########################################
########   class simulatedTS
#########################################

##' An S4 class to store a list of DAGs and associated observations
##' @name simulatedTS
##' @param list.DAGs : list of stored DAGs
##' @param list.observationsDAGs : list of observed datasets, each sampled from the corresponding member of list.DAGs
##' @param NDAG  : number of DAGs.
##' @param functionType : type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid")
##' @param seed : random seed
setClass("simulatedTS",
         slots = list(list.DAGs="list",list.observationsDAGs="list",
                      NDAG="numeric", list.numTS="list",
                      list.fs="list", list.Y="list", list.doNeigh="list",noNodes="numeric"))

##' creation of a "simulatedTS" containing a list of DAGs and associated time series observations
##' @param .Object : simulatedTS object
##' @param NDAG : number of DAGs to be created and simulated
#' @param noNodes  : number of Nodes of the DAGs. If it is a two-valued vector , the value of Nodes is randomly sampled in the interval
#' @param N  : number of sampled observations for each DAG. If it is a two-valued vector [a,b], the value of N is randomly sampled in the interval [a,b]
#' @param sdn : standard deviation of aditive noise. If it is a two-valued vector, the value of N is randomly sampled in the interval
#' @param typeser : type time series: numbers associated to different processes in function genSTAR 
#' @param seed : random seed
#' @param verbose : if TRUE it prints out the state of progress
#' @param goParallel : if TRUE it uses parallelism
#' @param nseries : number of series (if > 1 it is multivariate)
#' @references Gianluca Bontempi, Maxime Flauder (2015) From dependency to causality: a machine learning approach. JMLR, 2015, \url{http://jmlr.org/papers/v16/bontempi15a.html}
#' @examples
#' require(RBGL)
#' require(gRbase)
#' require(foreach)
#' descr=new("D2C.descriptor")
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedTS",NDAG=10, N=c(50,100),noNodes=c(15,40), typser=3,
#'              functionType = "linear", seed=0,sdn=c(0.45,0.75))
#' @export
#'
#'
setMethod("initialize",
          "simulatedTS",
          function(.Object, NDAG=1,
                   noNodes=sample(10:20,size=1),
                   verbose=TRUE,N=sample(100:500,size=1),
                   typeser=1:5,
                   seed=1234,sdn=0.5, goParallel=FALSE, nseries=1)
          {
            
            ##generate a training set
            ## NDAG the number of network to use
            
            
            if (goParallel){
              # cl <- makeForkCluster(10)
              #  registerDoParallel(cl)
              `%op%` <-  `%dopar%` 
            }   else {
              `%op%` <-`%do%`
            }
            
            
            X=NULL
            Y=NULL
            list.DAGs=NULL
            list.observationsDAGs=NULL
            if (NDAG<=0)
              return(.Object)
            
            FF<-foreach (i=1:NDAG) %op%{
              ##        for (i in 1:NDAG){
              set.seed(seed+i)
              
              nseries.i<-nseries
              if (length(nseries)>1)
                nseries.i<-sample(nseries[1]:nseries[2],1)
              
              N.i<-N
              if (length(N)>1)
                N.i<-sample(N[1]:N[2],1)
              
              noNodes.i<-max(4,noNodes[1])
              if (length(noNodes)==2)
                noNodes.i<-sample(max(4:noNodes[1]):max(4:noNodes[2]),1)
              
              
              sdn.i<-sdn
              if (length(sdn)>1)
                sdn.i<-runif(1,sdn[1],sdn[2])
              
              
              num=sample(typeser,1)
              
              if (nseries.i>1)
                G<-genSTAR(n=nseries.i,nn=noNodes.i,NN=N.i,sd=sdn.i,num=num,loc=1)
              else
                G<-genTS(nn=noNodes.i,N=N.i,sd=sdn.i,num=num)
              
              if (any(is.na(G$D) | is.infinite(G$D)))
                stop("error in data generation")
              netwDAG<-G$DAG 
              nodes(netwDAG)<-as.character(1:NCOL(G$D))
              observationsDAG = G$D
              fsTS=G$fs
              Y=G$Y
              doNeigh=G$doNeigh
              
              if (verbose){
                
                cat("simulatedTS: TS number:",i,":",num," generated: #nodes=", NCOL(observationsDAG),
                    "# edges=",sum(unlist(lapply(graph::edges(netwDAG),length))), "# samples=", N.i, "\n")
                
              }
              
              list(observationsDAG=observationsDAG,netwDAG=netwDAG,
                   numTS=num,fsTS=fsTS,Y=Y,doNeigh=doNeigh)
            } ## foreach
            
            
            .Object@list.DAGs=lapply(FF,"[[",2)
            .Object@list.observationsDAGs=lapply(FF,"[[",1)
            .Object@list.numTS=lapply(FF,"[[",3)
            .Object@list.fs=lapply(FF,"[[",4)
            .Object@list.Y=lapply(FF,"[[",5)
            .Object@list.doNeigh=lapply(FF,"[[",6)
            to.remove=which(unlist(lapply(lapply(.Object@list.DAGs,edgeList),length))==0)
            if (length(to.remove)>0){
              .Object@list.DAGs=.Object@list.DAGs[-to.remove]
              .Object@list.observationsDAGs=.Object@list.observationsDAGs[-to.remove]
              .Object@list.numTS=Object@list.numTS[-to.remove]
              .Object@list.fs=Object@list.fs[-to.remove]
              .Object@list.Y=Object@list.Y[-to.remove]
              Object@list.doNeigh=Object@list.doNeigh[-to.remove]
            }
            .Object@NDAG=length(.Object@list.DAGs)
            .Object@noNodes=noNodes
            .Object
          }
)



setGeneric("update", def=function(object,...) {standardGeneric("update")})

#' update of a "simulatedDAG" with a list of DAGs and associated observations
#' @param object :  simulatedDAG to be updated
#' @param list.DAGs : list of stored DAGs
#' @param list.observationsDAGs : list of observed datasets, each sampled from the corresponding member of list.DAGs
#' @export
setMethod(f="update",
          signature="simulatedDAG",
          definition=function(object,list.DAGs,list.observationsDAGs) {
            if (length(list.DAGs)!=length(list.observationsDAGs))
              stop("Lists with different lengths !")
            object@list.DAGs=c(object@list.DAGs,list.DAGs)
            object@list.observationsDAGs=c(object@list.observationsDAGs,list.observationsDAGs)
            object@NDAG=length(object@list.DAGs)
            object
            
          }
)



#########################################
########   class D2C
#########################################

setOldClass("randomForest")

#' An S4 class to store the RF model trained on the basis of the descriptors of NDAG DAGs
setClass("D2C",
         slots = list(#mod="randomForest", 
           mod="list", X="matrix",origX="matrix",Y="numeric",
           descr="D2C.descriptor",scaled="numeric",features="numeric",center="numeric",
           allEdges="list",ratioMissingNode="numeric",ratioEdges="numeric",
           max.features="numeric",type="character",classifier="character"
         ))

#' creation of a D2C object which preprocesses the list of DAGs and observations contained in sDAG and fits a  Random Forest classifier
#' @name  D2C object
#' @param .Object : the D2C object
#' @param sDAG : simulateDAG object
#' @param descr  : D2C.descriptor object containing the parameters of the descriptor
#' @param max.features  : maximum number of features used by the Random Forest classifier \link[randomForest]{randomForest}. The features are selected by the importance returned by the function \link[randomForest]{importance}.
#' @param ratioEdges  : percentage of existing edges which are added to the training set
#' @param ratioMissingNode  : percentage of existing nodes which are not considered. This is used to emulate latent variables.
#' @param goParallel : if TRUE it uses parallelism
#' @param verbose  : if TRUE it prints the state of progress
#' @param type : type of predicted dependency. It takes values in \{ \code{is.parent, is.child, is.ancestor, is.descendant, is.mb} \}
#' @param EErep: Easy Ensemble size to deal with unbalancedness
#' @param  rev: if TRUE, it uses both directions of the edge to train the learner (i.e. if i is parent of j, it is also true that j is a child of i)
#' @references Gianluca Bontempi, Maxime Flauder (2015) From dependency to causality: a machine learning approach. JMLR, 2015, \url{http://jmlr.org/papers/v16/bontempi15a.html}
#' @examples
#' require(RBGL)
#' require(gRbase)
#'  require(foreach)
#' descr=new("D2C.descriptor")
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedDAG",NDAG=2, N=50,noNodes=10,
#'              functionType = "linear", seed=0,sdn=0.5)
#' example<-new("D2C",sDAG=trainDAG, descr=descr.example)
#' @export
setMethod("initialize",
          "D2C",
          function(.Object, sDAG,
                   descr=new("D2C.descriptor"),
                   verbose=TRUE,
                   ratioMissingNode=0,
                   ratioEdges=1,max.features=20,
                   goParallel=FALSE,npar=5,
                   type="is.parent", rev=TRUE) {
            
            #generate a training set
            # NDAG the number of network to use
            #functionType example : "R1" "R2" "sigmoid1"
            `%op%` <- if (goParallel) `%dopar%` else `%do%`
            
            
            .Object@descr=descr
            .Object@ratioMissingNode= ratioMissingNode
            .Object@ratioEdges= ratioEdges
            .Object@max.features= max.features
            .Object@type= type
            
            X=NULL
            Y=NULL
            allEdges=NULL
            FF<-NULL
            
            
            iter=1
            while ( iter <=sDAG@NDAG){
              iFF<-foreach (ii=iter:min(sDAG@NDAG,iter+npar-1)) %dopar%{
                ##  FF<-foreach (ii=1:sDAG@NDAG) %op%{
                ## for (ii in iter:min(sDAG@NDAG,iter+npar-1))  {   ### D2C
                
                set.seed(ii)
                
                DAG = sDAG@list.DAGs[[ii]]
                observationsDAG =sDAG@list.observationsDAGs[[ii]]
                if (!is.vector(observationsDAG)){
                  if (NROW(observationsDAG)<30)
                    observationsDAG=observationsDAG[sample(NROW(observationsDAG),50,rep=TRUE),]
                  if (verbose)
                    cat("D2C:  DAG", ii, "/", sDAG@NDAG, "(N,n)=", dim(observationsDAG), " ")
                  if (any(apply(observationsDAG,2,sd)<0.001))
                    cat(" constant values in DAG data ")
                  
                  Nodes = nodes(DAG)
                  
                  sz=max(2,ceiling(length(Nodes)*(1-ratioMissingNode)))
                  keepNode = sort(sample(Nodes,
                                         size = sz ,
                                         replace = F))
                  
                  DAG2 =subGraph(keepNode, DAG) ## subGraph {graph}
                  
                  
                  iDAG2=graph.adjacency(as(DAG2,"matrix"))  ## as(DAG2,"matrix"): adjacency matrix of graphNEL DAG
                  ##  transforms the graphNEL adjacency matrix into igraph object
                  # Use as_graphnel to transform an igraph into graphNEL
                  
                  ##choose which edge to train / predict and find the right label
                  nEdge = length(edgeList(DAG))
                  sz=max(1,round(nEdge*ratioEdges))
                  N0=0
                  N1=0
                  cnt=0
                  while(cnt < 2 & (N0<10 | N1<10)){
                      cnt=cnt+1
                      edgesM = matrix(unlist(sample(edgeList(DAG2),
                                                    size = sz,replace = F)),ncol=2,byrow = TRUE)
                      edgesM = rbind(edgesM,t(replicate(n =2*sz ,
                                                        sample(keepNode,size=2,replace = FALSE)))) ## random edges
                     
                    if (type!="is.parent") {
                      edgesM = rbind(edgesM,t(replicate(n =2*sz ,
                                           sample(keepNode,size=2,replace = FALSE)))) ## random edges
                    }
                    nEdges =  NROW(edgesM)
                    
                    if (verbose)
                      cat("nEdges=", nEdges, " ")
                    
                    
                    
                    if (rev)
                      labelEdge = numeric(2*nEdges)
                    else
                      labelEdge = NULL
                    
                    ##compute the descriptor for the edges
                    X.out = NULL
                    
                    if (rev){
                      for(j in 1:nEdges){
                        
                        I =as(edgesM[j,1],"numeric") 
                        J =as(edgesM[j,2],"numeric") 
                        
                        d<-descriptor(observationsDAG,I,J,lin=descr@lin,acc=descr@acc,
                                      struct=descr@struct,bivariate=descr@bivariate,
                                      pq=descr@pq,ns=descr@ns,maxs=descr@maxs,boot=descr@boot,
                                      errd=descr@residual, delta=descr@diff, stabD=descr@stabD)
                        
                        
                        if (type=="is.parent")
                          if (is.parent(iDAG2,edgesM[j,1],edgesM[j,2]))
                            labelEdge[(2*j)-1] =1
                        if (type=="is.child")
                          if (is.child(iDAG2,edgesM[j,1],edgesM[j,2]))
                            labelEdge[(2*j)-1] =1
                        if (type=="is.ancestor")
                          if (is.ancestor(iDAG2,edgesM[j,1],edgesM[j,2]))
                            labelEdge[(2*j)-1] =1
                        if (type=="is.descendant")
                          if (is.descendant(iDAG2,edgesM[j,1],edgesM[j,2]))
                            labelEdge[(2*j)-1] =1
                        if (type=="is.mb")
                          if (is.mb(iDAG2,edgesM[j,1],edgesM[j,2]))
                            labelEdge[(2*j)-1] =1
                        X.out = rbind(X.out,d)
                        
                        ## reverse edge
                        I =as(edgesM[j,2],"numeric") ;
                        J =as(edgesM[j,1],"numeric") ;
                        
                        
                        d<-descriptor(observationsDAG,I,J,lin=descr@lin,acc=descr@acc,
                                      struct=descr@struct,bivariate=descr@bivariate,
                                      pq=descr@pq,ns=descr@ns,maxs=descr@maxs,boot=descr@boot,
                                      errd=descr@residual, delta=descr@diff,stabD=descr@stabD)
                        
                        
                        if (type=="is.parent")
                          if (is.parent(iDAG2,edgesM[j,2],edgesM[j,1]))
                            labelEdge[2*j] =1
                        if (type=="is.child")
                          if (is.child(iDAG2,edgesM[j,2],edgesM[j,1]))
                            labelEdge[2*j] =1
                        if (type=="is.ancestor")
                          if (is.ancestor(iDAG2,edgesM[j,2],edgesM[j,1]))
                            labelEdge[2*j] =1
                        if (type=="is.descendant")
                          if (is.descendant(iDAG2,edgesM[j,2],edgesM[j,1]))
                            labelEdge[2*j] =1
                        if (type=="is.mb")
                          if (is.mb(iDAG2,edgesM[j,2],edgesM[j,1]))
                            labelEdge[2*j] =1
                        X.out = rbind(X.out,d)
                        
                      }
                    } else {    ### if rev
                      for (j in 1:nEdges){
                        I =as(edgesM[j,1],"numeric") ; # parent
                        J =as(edgesM[j,2],"numeric") ; # child
                        fs<-timecauses(NCOL(observationsDAG),sDAG@noNodes,J)
                        if (is.element(I,fs)){
                          dfs<-unique(c(I,J,fs))
                          ## it considers the edge only if this is temporally feasible
                          d<-descriptor(observationsDAG,I,J,lin=descr@lin,acc=descr@acc,
                                        struct=descr@struct,bivariate=descr@bivariate,
                                        pq=descr@pq,ns=descr@ns,maxs=descr@maxs,boot=descr@boot,
                                        errd=descr@residual, delta=descr@diff,stabD=descr@stabD)
                          
                          
                          if (type=="is.parent")
                            if (is.parent(iDAG2,edgesM[j,1],edgesM[j,2])){
                              labelEdge =c(labelEdge,1)
                            }else{
                              labelEdge =c(labelEdge,0)
                            }
                          if (type=="is.child")
                            if (is.child(iDAG2,edgesM[j,1],edgesM[j,2])){
                              labelEdge =c(labelEdge,1)
                            }else{
                              labelEdge =c(labelEdge,0)
                            }
                          if (type=="is.ancestor")
                            if (is.ancestor(iDAG2,edgesM[j,1],edgesM[j,2])){
                              labelEdge =c(labelEdge,1)
                            }else{
                              labelEdge =c(labelEdge,0)
                            }
                          if (type=="is.descendant")
                            if (is.descendant(iDAG2,edgesM[j,1],edgesM[j,2])){
                              labelEdge =c(labelEdge,1)
                            }else{
                              labelEdge =c(labelEdge,0)
                            }
                          if (type=="is.mb")
                            if (is.mb(iDAG2,edgesM[j,1],edgesM[j,2])){
                              labelEdge =c(labelEdge,1)
                            }else{
                              labelEdge =c(labelEdge,0)
                            }
                          X.out = rbind(X.out,d)
                          
                        }
                      }  
                    } ## if rev
                    N0=length(which(labelEdge==0))
                    N1=length(which(labelEdge==1))
                  } ## while
                  if (verbose)
                    cat("Descriptor (N,n)=", dim(X.out), 
                        "N0=",length(which(labelEdge==0)), 
                        "N1=",length(which(labelEdge==1)),"DONE \n")
                  
                  list(X=X.out,Y=labelEdge,edges=edgesM)
                }  ## if (NROW(observationsDAG)>10)
              } ## foreach
              
              iter=iter+npar
              FF<-c(FF,iFF)
              if (verbose)
                cat(length(FF),"DAG processed \n")
              
            } ## while
            X<-do.call(rbind,lapply(FF,"[[",1))
            Y<-do.call(c,lapply(FF,"[[",2))
            allEdges<-lapply(FF,"[[",3)
            
            
            
            .Object@origX<-X
            .Object@Y=Y
            .Object@allEdges=allEdges
            return(.Object)
            
            
          }
)

#' @docType methods
setGeneric("makeModel", def=function(object,...) {standardGeneric("makeModel")})

#' creation of a D2C object which preprocesses the list of DAGs and observations contained in sDAG and fits a  Random Forest classifier
#' @name  D2C object
#' @param .Object : the D2C object
#' @param sDAG : simulateDAG object
#' @param descr  : D2C.descriptor object containing the parameters of the descriptor
#' @param max.features  : maximum number of features used by the Random Forest classifier \link[randomForest]{randomForest}. The features are selected by the importance returned by the function \link[randomForest]{importance}.
#' @param ratioEdges  : percentage of existing edges which are added to the training set
#' @param ratioMissingNode  : percentage of existing nodes which are not considered. This is used to emulate latent variables.
#' @param goParallel : if TRUE it uses parallelism
#' @param subset : if bivar it uses only bivariate descriptors
#' @param verbose  : if TRUE it prints the state of progress
#' @param EErep: Easy Ensemble size to deal with unbalancedness
#' @references Gianluca Bontempi, Maxime Flauder (2015) From dependency to causality: a machine learning approach. JMLR, 2015, \url{http://jmlr.org/papers/v16/bontempi15a.html}
#' @examples
#' @export
setMethod("makeModel",
          "D2C",
          function(object, max.features=30,
                   classifier="RF",
                   EErep=5,verbose=TRUE,subset="all") {
            
            X<-object@origX
            Y<-object@Y
            features=1:NCOL(X)
            
            if (subset=="bivar")
              features<- grep('B.',colnames(X))
            if (subset=="multivar")
              features<- grep('M.',colnames(X))
            if (subset=="noInt"){
              featInt<- grep('Int.',colnames(X))
              features<-setdiff(features,featInt)
            }
            
            X<-X[,features]
            wna<-which(apply(X,2,sd)<0.001)
            if (length(wna)>0){
              features<-setdiff(features,features[wna])
              X=X[,-wna] 
            }
            wna<-which(is.na(apply(X,2,mean)))
            if (length(wna)>0){
              features<-setdiff(features,features[wna])
              X=X[,-wna] 
            }
            
            X<-scale(X)
            object@scaled=attr(X,"scaled:scale")
            object@center=attr(X,"scaled:center")
            object@classifier=classifier
            
            object@features=features
            
            
            listRF<-list()
            #featrank<-mrmr(X ,factor(Y),min(NCOL(X),3*max.features))
            if (classifier=="RF"){
              RF <- randomForest(x =X ,y = factor(Y),importance=TRUE)
              IM<-importance(RF)[,"MeanDecreaseAccuracy"]
            }
            if (classifier=="XGB.1"){
              RF <- xgboost(data =X ,label = Y,nrounds=20,objective = "binary:logistic",eta=0.1)
              IM<-numeric(NCOL(X))
              names(IM)=colnames(X)
              IM[xgb.importance(model = RF)[,1]]=xgb.importance(model = RF)[,2]
            }
            if (classifier=="XGB.2"){
              RF <- xgboost(data =X ,label = Y,nrounds=20,objective = "binary:logistic",eta=0.2)
              IM<-numeric(NCOL(X))
              names(IM)=colnames(X)
              IM[xgb.importance(model = RF)[,1]]=xgb.importance(model = RF)[,2]
            }
            featrank<- sort(IM,decr=TRUE,ind=TRUE)$ix
            if (verbose)
              cat("Best descriptors: ", colnames(X)[featrank], "\n")
            
            ratio=3
            if (classifier=="RF")
              ratio=1
            
            
            for (rep in 1:EErep){
              w0<-which(Y==0)
              w1<-which(Y==1)
              if (length(w0)>=ratio*length(w1))
                w0<-sample(w0,round(ratio*length(w1)))
              
              if (length(w1)>=length(w0))
                w1<-sample(w1,round(length(w0)))
              Xb<-X[c(w0,w1),]
              Yb<-Y[c(w0,w1)]
              
              
              rank<-featrank
              rank<-rank[1:min(max.features,length(rank))]
              Xb=Xb[,rank]
              if (classifier=="RF")
                RF <- randomForest(x =Xb ,y = factor(Yb))
              if (classifier=="XGB.1")
                RF=xgboost(data =Xb ,label = Yb,nrounds=5,objective = "binary:logistic",eta=0.1)
              if (classifier=="XGB.2")
                RF=xgboost(data =Xb ,label = Yb,nrounds=5,objective = "binary:logistic",eta=0.2)
              
              listRF<-c(listRF,list(list(mod=RF,feat=rank)))
              if (verbose)
                cat(classifier, " ", rep, ":",RF$confusion, " : (N,n)=", dim(Xb), "\n")
            } ## for rep
            
            object@mod=listRF
            
            object
          }
)


#' predict if there is a connection between node i and node j
#' @param object : a D2C object
#' @param i :  index of putative cause (\eqn{1 \le i \le n})
#' @param j  : index of putative effect (\eqn{1 \le j \le n})
#' @param data : dataset of observations from the DAG
#' @return list with binary response and probability of the existence of a directed edge
#' @examples
#' require(RBGL)
#' require(gRbase)
#' require(foreach)
#' data(example)
#'## load the D2C object
#' testDAG<-new("simulatedDAG",NDAG=1, N=50,noNodes=5,
#'            functionType = "linear", seed=1,sdn=c(0.25,0.5))
#' ## creates a simulatedDAG object for testing
#' plot(testDAG@@list.DAGs[[1]])
#' ## plot the topology of the simulatedDAG
#' predict(example,1,2, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 1->2 exists
#' predict(example,4,3, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 4->3 exists
#' predict(example,4,1, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 4->1 exists
#' @references Gianluca Bontempi, Maxime Flauder (2015) From dependency to causality: a machine learning approach. JMLR, 2015, \url{http://jmlr.org/papers/v16/bontempi15a.html}
#' @export
setMethod("predict", signature="D2C",
          function(object,i,j,data, rep=1){ 
            out = list()
            
            if (any(apply(data,2,sd)<0.01))
              stop("Error in D2C::predict: Remove constant variables from dataset. ")
            Response<-NULL
            Prob<-NULL
            
            for (repS in 1:rep){
              ## repetition with different subsets of other variables
              others=setdiff(1:NCOL(data),c(i,j))
              if (repS>1)
                others=sample(others, round(2*length(others)/3))
              D=data[,c(i,j,others)]
              # move the concerned variables to the first two places
              
              #X_descriptor = descriptor(data,i,j,
              X_descriptor = descriptor(D,1,2,
                                        lin = object@descr@lin,
                                        acc = object@descr@acc,
                                        ns=object@descr@ns,
                                        maxs=object@descr@maxs,
                                        struct = object@descr@struct,
                                        pq = object@descr@pq, 
                                        bivariate =object@descr@bivariate, 
                                        boot=object@descr@boot,
                                        errd=object@descr@residual, delta=object@descr@diff,
                                        stabD=object@descr@stabD)
              
              if (any(is.infinite(X_descriptor)))
                stop("Error in D2C::predict: infinite value ")
              X_descriptor=X_descriptor[object@features]
              
              X_descriptor=scale(array(X_descriptor,c(1,length(X_descriptor))),
                                 object@center,object@scaled)
              if (any(is.infinite(X_descriptor)))
                stop("error in D2C::predict")
              
              for (r in 1:length(object@mod)){
                mod=object@mod[[r]]$mod
                fs=object@mod[[r]]$feat
                #Response = c( Response, predict(mod, X_descriptor[fs], type="response"))
                if (object@classifier=="RF")
                  Prob = c(Prob,predict(mod, X_descriptor[fs], type="prob")[,"1"])
                if (length(grep("XGB",object@classifier))>=1)
                  Prob = c(Prob,predict(mod, array(X_descriptor[fs],c(1,length(fs)))))
                
              }
            }
            out[["response"]] =round(mean(Prob))
            out[["prob"]]=mean(Prob)
            return(out)
          })


#' @docType methods
setGeneric("joinD2C", def=function(object,...) {standardGeneric("joinD2C")})

#' update of a "D2C" with a list of DAGs and associated observations
#' @name join current D2C and input D2C
#' @param object :  D2C to be updated
#' @param input :  D2C to be joined
#' @param verbose : TRUE or FALSE
#' @param goParallel : if TRUE it uses  parallelism
#' @export
setMethod(f="joinD2C",
          signature="D2C",
          definition=function(object,input,
                              verbose=TRUE, goParallel= FALSE){
            
            `%op%` <- if (goParallel) `%dopar%` else `%do%`
            
            X<-rbind(object@origX,input@origX)
            Y<-c(object@Y,input@Y)
            features<-intersect(object@features,input@features)
            
            
            features<-1:NCOL(X)
            wna<-which(apply(X,2,sd)<0.01)
            if (length(wna)>0)
              features<-setdiff(features,wna)
            object@origX=X
            X<-scale(X[,features])
            object@scaled=attr(X,"scaled:scale")
            object@center=attr(X,"scaled:center")
            
            object@features=features
            object@Y=Y
            max.features=object@max.features
            
            listRF<-list()
            for (rep in 1:10){
              w0<-which(Y==0)
              w1<-which(Y==1)
              if (length(w0)>length(w1))
                w0<-sample(w0,length(w1))
              
              if (length(w1)>length(w0))
                w1<-sample(w1,length(w0))
              Xb<-X[c(w0,w1),]
              Yb<-Y[c(w0,w1)]
              
              
              RF <- randomForest(x =Xb ,y = factor(Yb),importance=TRUE)
              IM<-importance(RF)[,"MeanDecreaseAccuracy"]
              rank<-sort(IM,decr=TRUE,ind=TRUE)$ix[1:min(max.features,NCOL(Xb))]
              Xb=Xb[,rank]
              RF <- randomForest(x =Xb ,y = factor(Yb))
              
              listRF<-c(listRF,list(list(mod=RF,feat=rank)))
            }
            
            object@mod=listRF
            
            object
          }
)




#' @docType methods
setGeneric("updateD2C", def=function(object,...) {standardGeneric("updateD2C")})

#' update of a "D2C" with a list of DAGs and associated observations
#' @name update D2C
#' @param object :  D2C to be updated
#' @param sDAG : simulatedDAG object to update D2C
#' @param verbose : TRUE or FALSE
#' @param goParallel : if TRUE it uses  parallelism
#' @export
setMethod(f="updateD2C",
          signature="D2C",
          definition=function(object,sDAG,
                              verbose=TRUE, goParallel= FALSE){
            
            `%op%` <- if (goParallel) `%dopar%` else `%do%`
            ratioMissingNode=object@ratioMissingNode
            ratioEdges=object@ratioEdges
            descr=object@descr
            FF<-foreach (i=1:sDAG@NDAG) %op%{
              
              set.seed(i)
              DAG = sDAG@list.DAGs[[i]]
              observationsDAG =sDAG@list.observationsDAGs[[i]]
              
              Nodes = nodes(DAG)
              
              sz=max(2,ceiling(length(Nodes)*(1-ratioMissingNode)))
              keepNode = sort(sample(Nodes,
                                     size = sz ,
                                     replace = F))
              
              DAG2 =subGraph(keepNode, DAG)
              
              
              ##choose wich edge to train / predict and find the right label
              nEdge = length(edgeList(DAG))
              sz=max(1,round(nEdge*ratioEdges))
              
              edgesM = matrix(unlist(sample(edgeList(DAG2),
                                            size = sz,replace = F)),ncol=2,byrow = TRUE)
              edgesM = rbind(edgesM,t(replicate(n =sz ,
                                                sample(keepNode,size=2,replace = FALSE))))
              
              nEdges =  NROW(edgesM)
              labelEdge = numeric(nEdges)
              for(j in 1:nEdges){
                I =edgesM[j,1] ;
                J =edgesM[j,2] ;
                labelEdge[j] = as.numeric(I %in% inEdges(node = J,DAG2)[[1]])
              }
              
              
              ##compute the descriptor for the edges
              nNodes = length(labelEdge)
              
              X.out = NULL
              for(j in 1:nNodes){
                I =as(edgesM[j,1],"numeric") ;
                J =as(edgesM[j,2],"numeric") ;
                
                
                d<-descriptor(observationsDAG,I,J,lin=descr@lin,acc=descr@acc,
                              struct=descr@struct,bivariate=descr@bivariate,
                              pq=descr@pq,maxs=descr@maxs,ns=descr@ns,boot=descr@boot,
                              errd=descr@residual, delta=descr@diff, stabD=descr@stabD)
                
                
                
                
                
                X.out = rbind(X.out,d)
              }
              if (verbose)
                cat("D2C:  DAG", i, " processed \n")
              
              list(X=X.out,Y=labelEdge,edges=edgesM)
              
            }
            
            X<-do.call(rbind,lapply(FF,"[[",1))
            Y<-do.call(c,lapply(FF,"[[",2))
            allEdges<-lapply(FF,"[[",3)
            
            
            
            X<-scale(X[,object@features],attr(object@X,"scaled:center"),attr(object@X,"scaled:scale"))
            
            object@X=rbind(object@X,X)
            object@Y=c(object@Y,Y)
            object@allEdges=c(object@allEdges,allEdges)
            RF <- randomForest(x =object@X ,y = factor(object@Y),importance=TRUE)
            IM<-importance(RF)[,"MeanDecreaseAccuracy"]
            rank<-sort(IM,decr=TRUE,ind=TRUE)$ix[1:min(object@max.features,NCOL(X))]
            RF <- randomForest(x =object@X[,rank] ,y = factor(object@Y))
            object@rank=rank
            object@mod=RF
            
            object
            
          }
)

#' Dataset example
#'@title stored D2C object
#'@description small D2C object for testing D2C functionalities
#' @name example
#' @docType data
#' @keywords data
#' @export
#' @examples
#' require(RBGL)
#' require(gRbase)
#' data(example)
#' print(example@@mod)
#' ## Random Forest
#' print(dim(example@@X))
#' ## dimension of the training set
NULL




#' Dataset alarm
#'@title Alarm dataset
#'@description contains the adjacency matrix of the Alarm DAG (\code{true.net}) and the related measured dataset (\code{dataset}). See the vignette for an utilization of the dataset
#' @name alarm
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL

#' Benchmark alarm
#'@title Alarm benchmark
#'@description contains the adjacency matrix of the Alarm DAG (\code{true.net}) and the related measured dataset (\code{dataset}). See the vignette for an utilization of the dataset
#' @name alarm
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL


#' Adjacency matrix of the Alarm benchmark
#'@title Adjacency matrix of the Alarm dataset
#'@description contains the adjacency matrix of the Alarm DAG. See the vignette for an utilization of the dataset
#' @name true.net
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL

#'  Dataset of the Alarm benchmark
#'@title Dataset of the Alarm benchmark
#'@description contains the  measured dataset. See the vignette for an utilization of the dataset
#' @name dataset
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL
