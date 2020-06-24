library(D2C)
library(Rgraphviz)
graphics.off()

genSTAR<-function(n, nn,NN,sd=0.5,num=1,loc=2,verbose=FALSE){
  # n: number of time series
  # nn: max lag
  ## NN: number of samples
  ## loc : size neghborhood
  ##
  ## Returns:
  ## D [NN,n*nn] observation dataset
  ## DAG associated DAG
  if (nn<4)
    stop("Too few lags")
  Y=array(rnorm(n*nn,sd=0.1),c(nn,n))
  ep=0
  eold=numeric(n)
  eold2=numeric(n)
  th0=rnorm(1)
  fs<-sample(0:(nn-2),3)
  state=0
  print(num)
  for (ii in 1:NN){
    N=NROW(Y)
    y=numeric(n)
    if (num==1){
      nfs=2
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=-0.4*(3-mean(Y[N-fs[1],neigh])^2)/(1+mean(Y[N-fs[1],neigh])^2)+
          0.6*(3-(mean(Y[N-fs[2],neigh])-0.5)^3)/(1+(mean(Y[N-fs[2],neigh])-0.5)^4)+sd*e
        
      }
      
    }
    if (num==2){
      nfs=2
      
      #Y=c(Y,(0.4-2*exp(-50*Y[N-fs[1]]^2))*Y[N-fs[1]]+(0.5-0.5*exp(-50*Y[N-fs[2]]^2))*Y[N-fs[2]]+sd*(e+th0*ep))
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=(0.4-2*exp(-50*mean(Y[N-fs[1],neigh])^2))*mean(Y[N-fs[1],neigh])+
          (0.5-0.5*exp(-50*mean(Y[N-fs[2],neigh])^2))*mean(Y[N-fs[2],neigh])+sd*e
        
      }
    }
    if (num==3){
      nfs=3
      #Y=c(Y,1.5 *sin(pi/2*Y[N-fs[1]])-sin(pi/2*Y[N-fs[2]])+sin(pi/2*Y[N-fs[3]])+sd*(e+th0*ep))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=1.5 *sin(pi/2*mean(Y[N-fs[1],neigh]))-
          sin(pi/2*mean(Y[N-fs[2],neigh]))+
          sin(pi/2*mean(Y[N-fs[3],neigh]))+sd*e
        
      }
    }
    
    if (num==4){
      nfs=2
      #  Y=c(Y,2*exp(-0.1*Y[N-fs[1]]^2)*Y[N-fs[1]]-exp(-0.1*Y[N-fs[2]]^2)*Y[N-fs[2]]+sd*(e+th0*ep))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=2*exp(-0.1*mean(Y[N-fs[1],neigh])^2)*mean(Y[N-fs[1],neigh])-
          exp(-0.1*mean(Y[N-fs[2],neigh])^2)*mean(Y[N-fs[2],neigh])+sd*e
        
      }
    }
    
    if (num==5){
      nfs=1
      #Y=c(Y,-2*Y[N-fs[1]]*max(0,sign(-Y[N-fs[1]]))+0.4*Y[N-fs[1]]*max(0,sign(Y[N-fs[1]]))+sd*rnorm(1))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=-2*mean(Y[N-fs[1],neigh])*max(0,sign(-mean(Y[N-fs[1],neigh])))+
          0.4*mean(Y[N-fs[1],neigh])*max(0,sign(mean(Y[N-fs[1],neigh])))+sd*e
        
      }
    }
    if (num==6){
      nfs=2
      #Y=c(Y,0.8*log(1+3*Y[N-fs[1]]^2)-0.6*log(1+3*Y[N-fs[2]]^2)+sd*rnorm(1))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=0.8*log(1+3*mean(Y[N-fs[1],neigh])^2)-0.6*log(1+3*mean(Y[N-fs[2],neigh])^2)+sd*e
        
      }
    }
    if (num==7){
      nfs=2
      #y[i]=(0.4-2*cos(40*mean(Y[N-5,neigh]))*exp(-30*mean(Y[N-5,neigh])^2))*mean(Y[N-5,neigh])+(0.5-0.5*exp(-50*mean(Y[N-9,neigh])^2))*mean(Y[N-9,neigh])
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=(0.4-2*cos(40*mean(Y[N-fs[1],neigh]))*exp(-30*mean(Y[N-fs[1],neigh])^2))*mean(Y[N-fs[1],neigh])+
          (0.5-0.5*exp(-50*mean(Y[N-fs[2],neigh])^2))*mean(Y[N-fs[2],neigh])+sd*e
        
        
      }
    }
    if (num==8){
      nfs=2
      #Y=c(Y,(0.5-1.1*exp(-50*Y[N-fs[1]]^2))*Y[N]+(0.3-0.5*exp(-50*Y[N-fs[2]]^2))*Y[N-fs[2]]+sd*rnorm(1))
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=0.5-1.1*exp(-50*mean(Y[N-fs[1],neigh])^2)*mean(Y[N-fs[1],neigh])+
          (0.3-0.5*exp(-50*mean(Y[N-fs[2],neigh])^2))*mean(Y[N-fs[2],neigh])+sd*e
        
      }
      
    }
    
    if (num==9){
      nfs=2
      ##Y=c(Y,0.3*Y[N-fs[1]]+0.6*Y[N-fs[2]]+(0.1-0.9*Y[N-fs[1]]+0.8*Y[N-fs[2]])/(1+exp(-10*Y[N-fs[1]]))+sd*rnorm(1))
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=0.3*mean(Y[N-fs[1],neigh])+0.6*mean(Y[N-fs[2],neigh])+
          (0.1-0.9*mean(Y[N-fs[1],neigh])+
             0.8*mean(Y[N-fs[2],neigh]))/(1+exp(-10*mean(Y[N-fs[1],neigh])))+sd*e
      }
      
    }
    if (num==10){
      nfs=1
      ##Y=c(Y,sign(Y[N-fs[1]])+sd*rnorm(1))
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=sign(mean(Y[N-fs[1],neigh]))+sd*e
        
      }
    }
    if (num==11){
      nfs=1
      ##Y=c(Y,0.8*Y[N-fs[1]]-0.8*Y[N-fs[1]]/(1+exp(-10*Y[N-fs[1]]))+sd*rnorm(1)) 
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=0.8*mean(Y[N-fs[1],neigh])-0.8*mean(Y[N-fs[1],neigh])/(1+exp(-10*mean(Y[N-fs[1],neigh])))+sd*e
        
      }
    }
    if (num==12){
      nfs=2
      ##  Y=c(Y,0.3*Y[N-fs[1]]+0.6*Y[N-fs[2]]+(0.1-0.9*Y[N-fs[1]]+0.8*Y[N-fs[2]])/(1+exp(-10*Y[N-fs[1]]))+sd*rnorm(1))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=0.3*mean(Y[N-fs[1],neigh])+0.6*mean(Y[N-fs[2],neigh])+
          (0.1-0.9*mean(Y[N-fs[1],neigh])+0.8*mean(Y[N-fs[2],neigh]))/(1+exp(-10*mean(Y[N-fs[1],neigh])))+sd*e
        
      }
    }
    if (num==13){
      nfs=1
      ##  Y=c(Y,min(1,max(0,3.8*Y[N-fs[1]]*(1-Y[N-fs[1]])+sd*rnorm(1))))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=min(1,max(0,3.8*mean(Y[N-fs[1],neigh])*(1-mean(Y[N-fs[1],neigh]))+sd*e))
        
      }
    }
    if (num==14){
      nfs=2
      ##   Y=c(Y,1-1.4*Y[N-fs[1]]*Y[N-fs[1]] + 0.3*(Y[N-fs[2]])+0.001*sd*rnorm(1))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]= 1-1.4*mean(Y[N-fs[1],neigh])*mean(Y[N-fs[1],neigh]) + 
          0.3*mean(Y[N-fs[2],neigh])+0.0001*sd*e
        
      }
      
    }
    if (num==15){
      nfs=1
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        if (mean(Y[N-fs[1],neigh])<1){
          y[i]=-0.5*mean(Y[N-fs[1],neigh])+sd*e
        } else {
          y[i]=0.4*mean(Y[N-fs[1],neigh])+sd*e
        }
      }
    }
    
    if (num==16){
      nfs=1
      ##   Y=c(Y,1-1.4*Y[N-fs[1]]*Y[N-fs[1]] + 0.3*(Y[N-fs[2]])+0.001*sd*rnorm(1))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        if (abs(mean(Y[N-fs[1],neigh]))<=1){
          y[i]=0.9*mean(Y[N-fs[1],neigh])+sd*e
        } else {
          y[i]=-0.3*mean(Y[N-fs[1],neigh])+sd*e
        }
      }
    } 
    if (num==17){
      nfs=1
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        if (state==1){
          y[i]=-0.5*mean(Y[N-fs[1],neigh])+sd*rnorm(1)
        } else {
          y[i]=0.4*mean(Y[N-fs[1],neigh])+sd*rnorm(1)
        }
        if (runif(1)>0.9)
          state=1-state
      }
    }
    
    if (num==18){
      
      nfs=4
      #GARCH 
      ##Y=c(Y,sqrt(0.000019+0.846*((Y[N-fs[1]])^2+0.3*(Y[N-fs[2]])^2+0.2*(Y[N-fs[3]])^2+0.1*(Y[N-fs[4]])^2) )*sd*rnorm(1))
      
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e=rnorm(1)
        y[i]=sqrt(0.000019+0.846*((mean(Y[N-fs[1],neigh]))^2+
                                    0.3*(mean(Y[N-fs[2],neigh]))^2+
                                    0.2*(mean(Y[N-fs[3],neigh]))^2+
                                    0.1*(mean(Y[N-fs[4],neigh]))^2) )*sd*e
        
        
      }
    }
    
    if (num==19){
      nfs=1
      e=numeric(n)
      # y=0.7 y(t-1)*e(t-2)+e(t)
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e[i]=rnorm(1)
        y[i]=0.7*(mean(Y[N-fs[1],neigh]))*sd*eold2[i]+sd*e[i]
      }
      eold2=eold
      eold=e
    }
    
    if (num==20){
      nfs=2
      e=numeric(n)
      # y=0.4 y(t-1)-0.3 y(t-2)+0.5*y(t-1)*e(t-1)+e(t)
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e[i]=rnorm(1)
        y[i]=0.4*(mean(Y[N-fs[1],neigh]))-0.3*(mean(Y[N-fs[2],neigh]))
        +0.5*(mean(Y[N-fs[1],neigh]))*sd*eold[i]+sd*e[i]
      }
      eold2=eold
      eold=e
    }
    
    if (num==21){
      nfs=1
      e=numeric(n)
      # y=0.7 abs(y(t-1))/(abs(y(t-1))+2)+e
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e[i]=rnorm(1)
        y[i]=0.7*abs(mean(Y[N-fs[1],neigh]))/ (abs(mean(Y[N-fs[1],neigh]))+2)+sd*e[i]
      }
      eold2=eold
      eold=e
    }
    
    if (num==22){
      nfs=1
      e=numeric(n)
      # y=0.9 y(t-1)- 0.8 y(t-1)/(1+exp(-10 y(t-1))) +e
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e[i]=rnorm(1)
        y[i]=0.9*mean(Y[N-fs[1],neigh])
        -0.8*mean(Y[N-fs[1],neigh])/ (1+exp(-10*mean(Y[N-fs[1],neigh])))+sd*e[i]
      }
      eold2=eold
      eold=e
    }
    if (num==23){
      nfs=2
      e=numeric(n)
      # y=0.3 y(t-1)+ 0.6 y(t-2)+ (0.1-0.9 y(t-1)+0.8 y(t-2))/(1+exp(-10 y(t-1))) +e
      for (i in 1:n){
        neigh=max(1,i-loc):min(n,i+loc)
        e[i]=rnorm(1)
        y[i]=0.3*mean(Y[N-fs[1],neigh]) +0.6 *mean(Y[N-fs[2],neigh])
        +(0.1-0.9*mean(Y[N-fs[1],neigh])+0.8*mean(Y[N-fs[2],neigh]))/ (1+exp(-10*mean(Y[N-fs[1],neigh])))+sd*e[i]
      }
      eold2=eold
      eold=e
    }
    
    Y<-rbind(Y,y)
    
    if (any(is.nan(Y) | abs(Y)>1000)){
      browser()
      stop("error")
      
    }
    
  } ## for i
  
  
  
  
  
  fs=fs[1:nfs]
  
  Y=scale(Y[nn:NROW(Y),])
  if (any(is.nan(Y) ))
    stop("error")
  #M=MakeEmbedded(Y,n=numeric(n)+nn,delay=numeric(n),hor=rep(1,n),w=1:n)
  netwDAG<-new("graphNEL", nodes=as.character(1:(n*nn)), edgemode="directed")
  
  fs=sort(fs)
  print(fs+1)
  for (i in 1:n){
    for (j in 1:nn){
      for (f in fs){
        if ((j+f)<nn){
          for (neigh in max(1,i-loc):min(n,i+loc)){
            netwDAG <- addEdge(as.character((neigh-1)*nn+j+f+1), as.character((i-1)*nn+j), netwDAG, 1)
            if (verbose)
              cat("i=",i,"j=",j,":",(neigh-1)*nn+j+f+1,"->",(i-1)*nn+j,"\n")
          }
        }
      }
    }
  }
  
  
  list(D=NULL,DAG=netwDAG,fs=fs)
}


G<-genSTAR(20, 4,100,sd=0.5,num=sample(setdiff(1:23,14)),loc=2,verb=TRUE)
print(G$fs+1)
print(dim(G$D))
print(G$DAG)
#plot(G$DAG)