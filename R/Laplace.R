Laplace<-function(adj){
  # dd = rowSums(adj)
  # ee = ifelse(dd>0,dd^(-0.5),0)
  # ee = diag(ee)
  # temp=ee%*%adj
  # L  = diag(1,n1,n1)-temp%*%ee
  p.temp=dim(adj)[1]
  dd = rowSums(adj)
  ee = ifelse(dd>0,dd^(-0.5),0)
  if(p.temp==1){
    ee=matrix(0,1,1)
  }else  ee = diag(ee)
  temp=ee%*%adj
  L= diag(1,(p.temp),(p.temp))-temp%*%ee

  return(L)
}