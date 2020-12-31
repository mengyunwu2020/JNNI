L_int<-function(adj.list){
  M=length(adj.list)
  p_temp=as.numeric(lapply(adj.list, function(x) dim(x)[1]))
  temp1=lapply(adj.list,adj.all)
  L=lapply(temp1,Laplace)
  p_hat=p_temp+p_temp*(p_temp-1)/2
  L1_int=matrix(0,(sum(p_hat)),(sum(p_hat)))
  s=1
  for (m in 1:M){
    L1_int[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]=L[[m]]
    s=s+p_hat[m]
  }

  return(L1_int)

}