library(pcalg)
p <- 5
rDAG <- randomDAG(p, prob = 0.2, lB=0.1, uB=1)

Admat<-gRbase::as.adjMAT(rDAG)

## generate 1000 samples of DAG using standard normal error distribution
n <- 100
d.normMat <- rmvDAG(n, rDAG, errDist="normal")

## generate 1000 samples of DAG using standard t(df=4) error distribution
d.t4Mat <- rmvDAG(n, rDAG, errDist="t4")

## generate 1000 samples of DAG using standard normal with a cauchy
## mixture of 30 percent
d.mixMat <- rmvDAG(n, rDAG, errDist="mix",mix=0.3)

